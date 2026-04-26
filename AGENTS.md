## Git policy

Do not create git commits.

You may inspect git status and diffs, for example:

- `git status`
- `git diff`
- `git diff --staged`
- `git log`

Do not run commands that create, amend, rewrite, push, or otherwise modify commits, including:

- `git commit`
- `git commit --amend`
- `git reset`
- `git rebase`
- `git merge`
- `git cherry-pick`
- `git tag`
- `git push`

Leave all changed files unstaged unless I explicitly ask you to stage them.
At the end of a task, summarize the changed files and suggest a commit message, but do not commit.
