# -*- coding: utf-8 -*-
# Script: generate_changelog.py
# date_created: 2025-06-04T12:11:19+02:00
# date_modified: 2025-06-04T12:11:23+02:00

"""
Description: Generate a changelog from git commit history
"""

import subprocess
from datetime import datetime

def get_git_log():
    """取得 git commit log"""
    result = subprocess.run(
        ["git", "log", "--pretty=format:%H%x09%an%x09%ad%x09%s", "--date=short"],
        stdout=subprocess.PIPE,
        text=True
    )
    return result.stdout.strip().split("\n")

def group_commits_by_date(commits):
    grouped = {}
    for line in commits:
        sha, author, date, message = line.split("\t")
        if date not in grouped:
            grouped[date] = []
        grouped[date].append((sha, author, message))
    return grouped

def generate_changelog(grouped):
    lines = ["# Changelog\n"]
    for date in sorted(grouped.keys(), reverse=True):
        lines.append(f"## {date}")
        for sha, author, message in grouped[date]:
            lines.append(f"- {message} ({author}, `{sha[:7]}`)")
        lines.append("")  # 空行
    return "\n".join(lines)

def write_changelog(content, filename="CHANGELOG.md"):
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)

if __name__ == "__main__":
    commits = get_git_log()
    grouped = group_commits_by_date(commits)
    changelog = generate_changelog(grouped)
    write_changelog(changelog)
    print("✅ CHANGELOG.md 已成功產生")
