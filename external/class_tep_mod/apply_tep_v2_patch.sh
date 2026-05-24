#!/bin/bash
#
# Apply TEP v2.0 Patch to CLASS
# Usage: ./apply_tep_v2_patch.sh [path_to_class]

set -e

# Get CLASS directory
CLASS_DIR=${1:-../class}
if [ ! -d "$CLASS_DIR" ]; then
    echo "Error: CLASS directory not found at $CLASS_DIR"
    echo "Usage: $0 /path/to/class_public"
    exit 1
fi

echo "========================================"
echo "TEP-CLASS v2.0 Patch Application"
echo "========================================"
echo ""
echo "Target CLASS directory: $CLASS_DIR"
echo ""

# Create backup
echo "Creating backup..."
BACKUP_DIR="${CLASS_DIR}_backup_$(date +%Y%m%d_%H%M%S)"
cp -r "$CLASS_DIR" "$BACKUP_DIR"
echo "Backup created: $BACKUP_DIR"
echo ""

# Check if already patched
grep -q "tep_mode" "$CLASS_DIR/include/background.h" 2>/dev/null && {
    echo "Warning: TEP parameters already found in background.h"
    echo "This directory may already be patched."
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
}

echo "========================================"
echo "Step 1: Modify background.h"
echo "========================================"

# Find the line with "};" at end of struct background
# We'll add TEP parameters before that line

BACKUP_H="$CLASS_DIR/include/background.h.backup"
cp "$CLASS_DIR/include/background.h" "$BACKUP_H"

# Create temporary file with TEP parameters
cat > /tmp/tep_params.txt << 'EOF'

  /* TEP (Temporal Equivalence Principle) parameters */
  short tep_mode;           /**< TEP mode flag (_TRUE_ or _FALSE_) */
  double epsilon_T;         /**< TEP coupling amplitude */
  double z_T;              /**< Characteristic redshift for TEP */
  double n_T;              /**< Power-law index for TEP transition */

EOF

echo "Adding TEP parameters to background structure..."

# Find line with "short has_*_density;" pattern (near end of struct background)
# Insert before the last few fields that are standard

# Use sed to add TEP parameters before the last field that ends with _scalar or similar marker
# This is a bit tricky, so we'll use a different approach

# Create a Python script to do the modification properly
cat > /tmp/patch_background_h.py << 'PYEOF'
import re
import sys

input_file = sys.argv[1]
output_file = input_file  # Overwrite in place

with open(input_file, 'r') as f:
    content = f.read()

# Find the struct background and add TEP parameters
# Look for pattern: struct background { ... };

tep_params = '''\n  /* TEP (Temporal Equivalence Principle) parameters */
  short tep_mode;           /**< TEP mode flag (_TRUE_ or _FALSE_) */
  double epsilon_T;         /**< TEP coupling amplitude */
  double z_T;              /**< Characteristic redshift for TEP */
  double n_T;              /**< Power-law index for TEP transition */
'''

# Find a good insertion point - look for "short" declaration near end of struct
# Pattern: look for "short some_field_name;" followed by eventually "};"
# We'll add before the last few standard fields

# Find the last occurrence of "short has_" in the struct
lines = content.split('\n')
insert_idx = None
for i, line in enumerate(lines):
    if 'short has_' in line and '_density;' in line:
        insert_idx = i

if insert_idx is None:
    # Fallback: find last field before closing brace
    for i, line in enumerate(lines):
        if line.strip() == '};' and i > 1000:  # Make sure it's the main struct close
            insert_idx = i
            break

if insert_idx is None:
    print("Error: Could not find insertion point in background.h")
    sys.exit(1)

# Insert TEP parameters before the insert point
lines.insert(insert_idx, tep_params)

with open(output_file, 'w') as f:
    f.write('\n'.join(lines))

print("Successfully added TEP parameters to background.h")
PYEOF

python3 /tmp/patch_background_h.py "$CLASS_DIR/include/background.h"

echo "✓ background.h modified"
echo ""

echo "========================================"
echo "Step 2: Modify common.h"
echo "========================================"

# Backup and add TEP string definitions
BACKUP_COMMON="$CLASS_DIR/include/common.h.backup"
cp "$CLASS_DIR/include/common.h" "$BACKUP_COMMON"

# Find line with "#define _TRUE_ 1" and add after it
grep -n "#define _TRUE_ 1" "$CLASS_DIR/include/common.h" | tail -1 | cut -d: -f1 > /tmp/true_line.txt
TRUE_LINE=$(cat /tmp/true_line.txt)

if [ -n "$TRUE_LINE" ]; then
    # Add TEP defines after _TRUE_ definition
    sed -i.tmp "${TRUE_LINE}a\\
\\
/* TEP parameter strings */\\
#define TEP_MODE_STRING \"tep_mode\"\\
#define TEP_EPSILON_T_STRING \"tep_epsilon_T\"\\
#define TEP_Z_T_STRING \"tep_z_T\"\\
#define TEP_N_T_STRING \"tep_n_T\"" "$CLASS_DIR/include/common.h"
    rm -f "$CLASS_DIR/include/common.h.tmp"
    echo "✓ common.h modified"
else
    echo "Warning: Could not find _TRUE_ definition in common.h"
fi
echo ""

echo "========================================"
echo "Step 3: Modify input.c"
echo "========================================"

BACKUP_INPUT="$CLASS_DIR/source/input.c.backup"
cp "$CLASS_DIR/source/input.c" "$BACKUP_INPUT"

# Create Python script to patch input.c
cat > /tmp/patch_input_c.py << 'PYEOF'
import re
import sys

input_file = sys.argv[1]

with open(input_file, 'r') as f:
    content = f.read()

# 1. Add TEP parameter reading in input_read_parameters function
# Find pattern: "Omega0_dcdmdr = pba->Omega0_dcdmdr;"
# Add TEP reading after background parameters section

tep_reading = '''
  /* TEP parameters */
  class_read_string(TEP_MODE_STRING,string1);
  if (strcmp(string1,"yes") == 0) {
    pba->tep_mode = _TRUE_;
  }
  else {
    pba->tep_mode = _FALSE_;
  }

  if (pba->tep_mode == _TRUE_) {
    class_read_double(TEP_EPSILON_T_STRING,pba->epsilon_T);
    class_read_double(TEP_Z_T_STRING,pba->z_T);
    class_read_double(TEP_N_T_STRING,pba->n_T);
    
    class_test(pba->epsilon_T < 0, errmsg, "epsilon_T must be non-negative");
    class_test(pba->z_T <= 0, errmsg, "z_T must be positive");
    class_test(pba->n_T <= 0, errmsg, "n_T must be positive");
  }
  else {
    pba->epsilon_T = 0.0;
    pba->z_T = 1.0;
    pba->n_T = 1.0;
  }

'''

# Find a good insertion point - after Omega0_dcdmdr reading
pattern = r'(Omega0_dcdmdr = pba->Omega0_dcdmdr;)'
match = re.search(pattern, content)
if match:
    insert_pos = match.end()
    content = content[:insert_pos] + '\n' + tep_reading + content[insert_pos:]
    print("✓ Added TEP parameter reading")
else:
    print("Warning: Could not find insertion point for TEP reading")

# 2. Add default values in input_default_params function
# Find "pba->Omega0_dcdmdr = 1e-5;" or similar default

tep_defaults = '''
  pba->tep_mode = _FALSE_;
  pba->epsilon_T = 0.0;
  pba->z_T = 1.0;
  pba->n_T = 1.0;
'''

pattern2 = r'(pba->Omega0_dcdmdr = [^;]+;)'
match2 = re.search(pattern2, content)
if match2:
    insert_pos2 = match2.end()
    content = content[:insert_pos2] + '\n' + tep_defaults + content[insert_pos2:]
    print("✓ Added TEP default values")
else:
    print("Warning: Could not find insertion point for TEP defaults")

with open(input_file, 'w') as f:
    f.write(content)

print("Successfully modified input.c")
PYEOF

python3 /tmp/patch_input_c.py "$CLASS_DIR/source/input.c"
echo "✓ input.c modified"
echo ""

echo "========================================"
echo "Step 4: Modify background.c"
echo "========================================"

BACKUP_BG="$CLASS_DIR/source/background.c.backup"
cp "$CLASS_DIR/source/background.c" "$BACKUP_BG"

# Create Python script to patch background.c
cat > /tmp/patch_background_c.py << 'PYEOF'
import re
import sys

input_file = sys.argv[1]

with open(input_file, 'r') as f:
    content = f.read()

# 1. Add TEP helper functions at the beginning of the file (after includes)
tep_functions = '''
/* TEP helper functions */
double tep_f_transition(struct background *pba, double z) {
    if (pba->tep_mode == _FALSE_ || pba->epsilon_T == 0.0) {
        return 0.0;
    }
    
    double z_effective = z;
    if (pba->z_T > 0 && z > pba->z_T * 3.0) {
        z_effective = pba->z_T * 3.0;
    }
    
    double ratio = z_effective / pba->z_T;
    double exponent = pow(ratio, pba->n_T);
    
    return 1.0 - exp(-exponent);
}

double tep_gamma_factor(struct background *pba, double z) {
    if (pba->tep_mode == _FALSE_ || pba->epsilon_T == 0.0) {
        return 1.0;
    }
    double f_T = tep_f_transition(pba, z);
    return 1.0 + pba->epsilon_T * f_T;
}

'''

# Find first function definition and insert before it
pattern = r'(/\*.*\*/\s*\n\s*int background_init)'
match = re.search(pattern, content, re.DOTALL)
if match:
    insert_pos = match.start()
    content = content[:insert_pos] + tep_functions + content[insert_pos:]
    print("✓ Added TEP helper functions")
else:
    print("Warning: Could not find insertion point for TEP functions")

# 2. Modify background_functions to apply TEP
# Find where pvecback[pba->index_bg_H] is set
# Add TEP modification after it

tep_modification = '''
  /* TEP modification if enabled */
  if (pba->tep_mode == _TRUE_ && pba->epsilon_T != 0.0) {
    double z = 1.0/a - 1.0;
    double H_original = pvecback[pba->index_bg_H];
    double gamma = tep_gamma_factor(pba, z);
    pvecback[pba->index_bg_H] = H_original * gamma;
    
    /* Recompute conformal Hubble rate */
    pvecback[pba->index_bg_H * pba->index_bg_a] = a * pvecback[pba->index_bg_H];
  }
'''

# Find the Hubble rate assignment and add TEP modification after
# Look for pattern where H is computed from rho_tot
pattern_h = r'(pvecback\[pba->index_bg_H\] = sqrt\(rho_tot - pvecback[^;]+;\))'
match_h = re.search(pattern_h, content)
if match_h:
    insert_pos_h = match_h.end()
    content = content[:insert_pos_h] + '\n' + tep_modification + content[insert_pos_h:]
    print("✓ Added TEP Hubble modification")
else:
    print("Warning: Could not find Hubble rate assignment for TEP modification")

with open(input_file, 'w') as f:
    f.write(content)

print("Successfully modified background.c")
PYEOF

python3 /tmp/patch_background_c.py "$CLASS_DIR/source/background.c"
echo "✓ background.c modified"
echo ""

echo "========================================"
echo "Step 5: Create test parameter file"
echo "========================================"

cat > "$CLASS_DIR/test/tep_v2.ini" << 'EOF'
# TEP v2.0 Test Parameter File
# Uses best-fit parameters from Pantheon+ SNe analysis

root = output/tep_v2_

# TEP parameters
tep_mode = yes
tep_epsilon_T = 0.1742
tep_z_T = 5.0
tep_n_T = 1.0

# Background parameters (TEP best-fit)
H0 = 72.82
omega_b = 0.0224
omega_cdm = 0.2386
omega_ncdm = 1.989e-4  # Single massive neutrino

# Neutrino parameters
N_ur = 2.0328
N_ncdm = 1
m_ncdm = 0.06

# Curvature (flat)
Omega_k = 0.0

# Output
output = tCl,pCl,lCl,mPk
lensing = yes
l_max_scalars = 2500
l_max_tensors = 500

# Pivot scale
k_pivot = 0.05

# Precision
P_k_max_h/Mpc = 10.
z_max_pk = 10.
z_pk = 0., 1., 10., 1090.

# Age of universe
age = 13.8
EOF

echo "✓ Created test/tep_v2.ini"
echo ""

echo "========================================"
echo "Step 6: Create LCDM comparison file"
echo "========================================"

cat > "$CLASS_DIR/test/tep_v2_lcdm.ini" << 'EOF'
# LCDM comparison for TEP v2.0
# Same parameters but with tep_mode = no

root = output/tep_v2_lcdm_

# TEP disabled
tep_mode = no
tep_epsilon_T = 0.0
tep_z_T = 1.0
tep_n_T = 1.0

# Same background parameters as TEP
H0 = 72.82
omega_b = 0.0224
omega_cdm = 0.2386

# ... rest same as tep_v2.ini
N_ur = 2.0328
N_ncdm = 1
m_ncdm = 0.06
Omega_k = 0.0

output = tCl,pCl,lCl,mPk
lensing = yes
l_max_scalars = 2500
l_max_tensors = 500
k_pivot = 0.05
P_k_max_h/Mpc = 10.
z_max_pk = 10.
z_pk = 0., 1., 10., 1090.
age = 13.8
EOF

echo "✓ Created test/tep_v2_lcdm.ini"
echo ""

echo "========================================"
echo "Patch Application Complete!"
echo "========================================"
echo ""
echo "Next steps:"
echo "  1. cd $CLASS_DIR"
echo "  2. make clean"
echo "  3. make"
echo "  4. ./class test/tep_v2.ini"
echo "  5. ./class test/tep_v2_lcdm.ini"
echo ""
echo "If compilation fails, restore from backup:"
echo "  cp -r $BACKUP_DIR $CLASS_DIR"
echo ""
