# GitHub first steps

This repository is being prepared privately while the project is in progress.

## Five useful terms

| Term | Meaning in this project |
| --- | --- |
| Repository | The project folder, together with its saved history |
| README | The front page explaining the research question and current progress |
| Commit | A saved snapshot of selected changes, with a short message explaining them |
| Push | Uploading your local commits to the GitHub copy |
| Private | Access is restricted to you and people you explicitly grant access to |

[GitHub's repository guide](https://docs.github.com/en/repositories/creating-and-managing-repositories/about-repositories) explains these concepts further.

## Your first visit

Open the repository on GitHub while signed in. Read the README, open the glaucoma script under `scripts/`, and explore the documentation under `docs/`. The initial commit records the first portfolio snapshot. The GitHub copy remains private until its visibility is deliberately changed.

## Updating the project later

Research still happens in R on your computer. Git records files after you save them; it does not capture unsaved R objects or automatically upload your workspace.

This repository holds a portable copy of the working precheck script. If you continue editing a different copy on an external drive, those edits do not automatically appear here. Bring the reviewed changes into `scripts/` while retaining its portable path setup, and update the precheck record to reflect what you actually ran.

Once a file has been edited and saved, open a terminal in this repository. These read-only commands show what changed:

```sh
git status
git diff
```

For example, after editing only the README, you can save and upload that change with:

```sh
git add README.md
git commit -m "Clarify project overview"
git push
```

`git add` selects a file for the next snapshot. `git commit` saves it locally. `git push` sends that snapshot to GitHub. For future script or documentation updates, select those specific files instead.

Describe progress precisely: adding plotting code and reviewing its plot are different milestones. Update the status after reviewing the actual output.
