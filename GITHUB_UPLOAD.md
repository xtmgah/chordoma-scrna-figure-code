# GitHub Upload Instructions

The local Git repository is ready at:

```bash
/Users/zhangt8/Terminal/Codex/Chrodoma_scRNA/chordoma-scrna-figure-code
```

Suggested GitHub repository:

```text
https://github.com/xtmgah/chordoma-scrna-figure-code
```

## Create the Empty Repository

1. Go to GitHub and create a new repository under `xtmgah`.
2. Use repository name `chordoma-scrna-figure-code`.
3. Set visibility to public if the manuscript will cite the repository directly.
4. Do not initialize with a README, `.gitignore`, or license, because this local repository already includes the needed files.

## Push From This Computer

After the empty GitHub repository exists and GitHub authentication is available locally:

```bash
cd /Users/zhangt8/Terminal/Codex/Chrodoma_scRNA/chordoma-scrna-figure-code
git push -u origin main
```

If HTTPS authentication is not configured, either install/use GitHub CLI authentication:

```bash
gh auth login
git push -u origin main
```

or switch the remote to SSH after adding an SSH key to GitHub:

```bash
git remote set-url origin git@github.com:xtmgah/chordoma-scrna-figure-code.git
git push -u origin main
```
