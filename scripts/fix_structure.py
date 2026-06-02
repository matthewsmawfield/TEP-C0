import os
import json
import re

base_dir = "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/site"
comp_dir = os.path.join(base_dir, "components")

# Replacements to apply
# For 5_micro_macro.html:
def fix_5_micro(content):
    content = content.replace("<h2>5. The Micro-Macro Handshake</h2>", "<h1>5. The Micro-Macro Handshake</h1>")
    content = content.replace("<h3>5.", "<h2>5.")
    content = content.replace("</h3>", "</h2>") # this will replace all h3 tags, which is fine since they are section headers
    return content

def fix_6_pioneer(content):
    content = content.replace("<h2>6. Empirical Tests: Pioneer and Planck</h2>", "<h1>6. Empirical Tests: Pioneer and Planck</h1>")
    content = content.replace("<h3>6.", "<h2>6.")
    content = content.replace("</h3>", "</h2>")
    return content

def fix_7_disc(content):
    content = content.replace("<h1>5. ", "<h1>7. ")
    content = content.replace("<h2>5.", "<h2>7.")
    return content

def fix_8_conc(content):
    content = content.replace("<h1>6. ", "<h1>8. ")
    return content

def fix_9_ref(content):
    content = content.replace("<h1>7. ", "<h1>9. ")
    return content

def fix_10_rep(content):
    content = content.replace("<h1>8. ", "<h1>10. ")
    return content

# Read and modify files
files_to_process = [
    ("5_micro_macro.html", "5_micro_macro.html", fix_5_micro),
    ("6_pioneer_planck.html", "6_pioneer_planck.html", fix_6_pioneer),
    ("5_discussion.html", "7_discussion.html", fix_7_disc),
    ("6_conclusion.html", "8_conclusion.html", fix_8_conc),
    ("7_references.html", "9_references.html", fix_9_ref),
    ("8_reproducibility.html", "10_reproducibility.html", fix_10_rep)
]

for old_name, new_name, func in files_to_process:
    old_path = os.path.join(comp_dir, old_name)
    if os.path.exists(old_path):
        with open(old_path, "r") as f:
            content = f.read()
        
        content = func(content)
        
        new_path = os.path.join(comp_dir, new_name)
        with open(new_path, "w") as f:
            f.write(content)
        
        if old_name != new_name:
            os.remove(old_path)

# Update manifest.json
manifest_path = os.path.join(base_dir, "manifest.json")
with open(manifest_path, "r") as f:
    manifest = json.load(f)

for section in manifest["sections"]:
    if section["id"] == "section-discussion":
        section["file"] = "7_discussion.html"
    elif section["id"] == "section-conclusion":
        section["file"] = "8_conclusion.html"
    elif section["id"] == "section-references":
        section["file"] = "9_references.html"
    elif section["id"] == "section-reproducibility":
        section["file"] = "10_reproducibility.html"

with open(manifest_path, "w") as f:
    json.dump(manifest, f, indent=2)

print("Formatting fixed!")
