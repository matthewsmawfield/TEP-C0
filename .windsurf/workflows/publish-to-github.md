---
description: Publish release to GitHub with version tag
---

# Publish to GitHub Workflow

Automated workflow to commit all changes, push to master branch, and create a versioned release tag using the commit message from `version.txt`.

## Prerequisites

- All publication preparation steps completed (/prepare-publication workflow)
- `version.txt` contains the commit message (0-6 high-level bullets, latest changes only)
- Git remote configured for GitHub
- Working directory is clean (no uncommitted changes that shouldn't be published)

## Phase 1: Verify Ready to Publish

```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

echo "=== Checking repository status ==="
echo "Current branch:"
git branch --show-current

echo ""
echo "Modified files (should only be publication-related):"
git status --short

echo ""
echo "Version to be released:"
head -1 version.txt

echo ""
echo "Commit message preview (from version.txt):"
head -1 version.txt
echo ""
tail -n +2 version.txt
```

## Phase 2: Extract Version and Commit Message

```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

# Extract version tag from version.txt (e.g., "v0.5 (Tortola)" -> "v0.5")
VERSION_TAG=$(head -1 version.txt | grep -oE '^v[0-9]+\.[0-9]+' || echo "")

# Full commit message from version.txt
COMMIT_MSG_FILE=$(mktemp)
tail -n +1 version.txt > "$COMMIT_MSG_FILE"

echo "Version tag will be: $VERSION_TAG"
echo "Commit message:"
cat "$COMMIT_MSG_FILE"
```

## Phase 3: Stage All Changes

// turbo
```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

echo "=== Staging all changes ==="
git add -A

echo "Staged files:"
git diff --cached --stat
```

## Phase 4: Commit with version.txt Message

// turbo
```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

# Commit using version.txt as the message
git commit -F version.txt

echo "=== Commit successful ==="
git log -1 --oneline
```

## Phase 5: Push to Master Branch

// turbo
```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

echo "=== Pushing to origin master ==="
git push origin master

echo "Push complete."
```

## Phase 6: Create Version Tag

// turbo
```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

# Extract version tag (e.g., "v0.5")
VERSION_TAG=$(head -1 version.txt | grep -oE '^v[0-9]+\.[0-9]+')
CODENAME=$(head -1 version.txt | grep -oE '\([^)]+\)' | tr -d '()' || echo "Release")

echo "Creating annotated tag: $VERSION_TAG"
git tag -a "$VERSION_TAG" -m "${VERSION_TAG} (${CODENAME}) - Release"

echo "Tag created:"
git tag -l -n1 "$VERSION_TAG"
```

## Phase 7: Push Tag to GitHub

// turbo
```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

echo "=== Pushing tag to origin ==="
git push origin "$VERSION_TAG"

echo ""
echo "Tag push complete!"
echo "GitHub release URL: https://github.com/matthewsmawfield/TEP-WB/releases/tag/$VERSION_TAG"
```

## Phase 8: Verify Publication

```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

echo "=== Verification ==="
echo "Latest commit:"
git log -1 --oneline

echo ""
echo "Tags on GitHub:"
git ls-remote --tags origin | grep -E 'refs/tags/v[0-9]\.[0-9]+$' | tail -5

echo ""
echo "Repository status (should be clean):"
git status

echo ""
echo "========================================"
echo "PUBLICATION COMPLETE"
echo "========================================"
echo "Version: $(head -1 version.txt)"
echo "Commit: $(git rev-parse --short HEAD)"
echo "Tag: $(git describe --tags --exact-match 2>/dev/null || echo 'N/A')"
echo "GitHub: https://github.com/matthewsmawfield/TEP-WB"
echo "========================================"
```

## Alternative: One-Command Publish

For subsequent releases after initial setup:

// turbo
```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB" && \
VERSION_TAG=$(head -1 version.txt | grep -oE '^v[0-9]+\.[0-9]+') && \
git add -A && \
git commit -F version.txt && \
git push origin master && \
git tag -a "$VERSION_TAG" -m "$(head -1 version.txt)" && \
git push origin "$VERSION_TAG" && \
echo "Published $(head -1 version.txt)"
```

## Troubleshooting

### Uncommitted Changes Exist
```bash
# Check what's not staged
git status

# If files should not be committed, stash them
git stash push -m "WIP before publication"

# Or reset if they shouldn't be kept
git checkout -- <file>
```

### Tag Already Exists
```bash
# Check existing tags
git tag -l

# If tag exists and needs updating (force - use with caution)
VERSION_TAG=$(head -1 version.txt | grep -oE '^v[0-9]+\.[0-9]+')
git tag -d "$VERSION_TAG"
git push origin ":refs/tags/$VERSION_TAG"
git tag -a "$VERSION_TAG" -m "$(head -1 version.txt)"
git push origin "$VERSION_TAG"
```

### Push Rejected
```bash
# Pull latest changes first
git pull origin master --rebase

# Then retry push
git push origin master
```

### Commit Message Issues
```bash
# If version.txt format is wrong, check it:
head -5 version.txt

# Should be:
# vX.Y (Codename) - Brief summary
#
# - Change 1
# - Change 2
```

## Post-Publication Checklist

After successful push:

- [ ] GitHub shows latest commit at https://github.com/matthewsmawfield/TEP-WB
- [ ] Tag appears in releases: https://github.com/matthewsmawfield/TEP-WB/releases
- [ ] Tag message matches version.txt first line
- [ ] PDF is downloadable from GitHub release assets (if manually uploaded)
- [ ] Website reflects new version (if auto-deployed via GitHub Pages)

## Manual GitHub Release Notes (Optional)

After pushing the tag, create rich release notes on GitHub:

1. Go to https://github.com/matthewsmawfield/TEP-WB/releases
2. Click "Draft a new release"
3. Select the tag created by this workflow
4. Add release title: `vX.Y (Codename)`
5. Copy changelog from `version.txt` (0-6 high-level bullets, latest changes only)
6. Attach PDF: `13-TEP-WB-vX.Y-Codename.pdf`
7. Mark as pre-release if applicable
8. Publish release

## Files Referenced

- `version.txt` - Source of commit message and version tag
- GitHub repository: https://github.com/matthewsmawfield/TEP-WB
- Default branch: master
