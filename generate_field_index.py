#!/usr/bin/env python3
"""Regenerate the field table in damage-types.md and fields/index.md from fields/*.md frontmatter.

Usage:
    python3 docs/generate_field_index.py

Run this after adding, removing, or editing any file in docs/fields/.
The table in damage-types.md between <!-- BEGIN FIELD TABLE --> and <!-- END FIELD TABLE -->
is replaced in-place. fields/index.md is fully rewritten.
"""
import os, re, sys

DOCS_DIR  = os.path.dirname(os.path.abspath(__file__))
FIELDS_DIR = os.path.join(DOCS_DIR, 'fields')
DAMAGE_MD  = os.path.join(DOCS_DIR, 'damage-types.md')
INDEX_MD   = os.path.join(FIELDS_DIR, 'index.md')

TIER_ORDER = {'summary': 0, 'standard': 1, 'diagnostics': 2, 'deprecated': 3}

BEGIN = '<!-- BEGIN FIELD TABLE -->'
END   = '<!-- END FIELD TABLE -->'


def parse_frontmatter(path):
    """Return dict of YAML frontmatter key: value (string values only)."""
    text = open(path).read()
    m = re.match(r'^---\n(.*?)\n---\n', text, re.DOTALL)
    if not m:
        return {}
    fm = {}
    for line in m.group(1).splitlines():
        if ':' in line:
            k, _, v = line.partition(':')
            fm[k.strip()] = v.strip()
    return fm


def load_fields():
    fields = []
    for fname in sorted(os.listdir(FIELDS_DIR)):
        if fname == 'index.md' or not fname.endswith('.md'):
            continue
        fm = parse_frontmatter(os.path.join(FIELDS_DIR, fname))
        if not fm:
            print(f'  warning: no frontmatter in {fname}', file=sys.stderr)
            continue
        fields.append({
            'file':      fname,
            'title':     fm.get('title', fname.replace('.md', '')),
            'tier':      fm.get('tier', 'standard'),
            'estimand':  fm.get('estimand', ''),
            'stability': fm.get('stability', 'stable'),
        })
    fields.sort(key=lambda f: (TIER_ORDER.get(f['tier'], 1), f['title']))
    return fields


def make_damage_table(fields):
    rows = ['| Block | File | Tier | Description |', '|-------|------|------|-------------|']
    for f in fields:
        title = f'`{f["title"]}`' if '_' in f['title'] or f['title'][0].islower() else f['title']
        rows.append(
            f'| {title} | [{f["file"]}](./fields/{f["file"]}) | {f["tier"]} | {f["estimand"]} |'
        )
    return '\n'.join(rows)


def make_index_md(fields):
    rows = ['| File | Title | Tier | Estimand |', '|------|-------|------|----------|']
    for f in fields:
        rows.append(
            f'| [{f["file"]}](./{f["file"]}) | `{f["title"]}` | `{f["tier"]}` | {f["estimand"]} |'
        )
    table = '\n'.join(rows)
    return f"""---
type: libtaph-field-index
title: libtaph JSON field reference
---

# libtaph JSON field reference

Each file documents one JSON block emitted by `fqdup profile --json`.
YAML frontmatter encodes `tier`, `estimand`, `stability` — consumed by `fqdup profile --json-schema`.

{table}

## Tier meanings

| Tier | Emitted by default | Use |
|------|--------------------|-----|
| `summary` | yes | automated pipelines, quick status |
| `standard` | yes | full analysis |
| `diagnostics` | `--json-level full` only | debugging, method development |
| `deprecated` | yes, with warning | legacy; removed in v3 |
"""


def update_damage_md(table_str):
    content = open(DAMAGE_MD).read()
    if BEGIN not in content or END not in content:
        sys.exit(f'error: sentinels not found in {DAMAGE_MD}')
    before = content[:content.index(BEGIN)]
    after  = content[content.index(END) + len(END):]
    open(DAMAGE_MD, 'w').write(before + BEGIN + '\n' + table_str + '\n' + END + after)


def main():
    fields = load_fields()
    print(f'Found {len(fields)} field files')

    table = make_damage_table(fields)
    update_damage_md(table)
    print(f'Updated {DAMAGE_MD}')

    open(INDEX_MD, 'w').write(make_index_md(fields))
    print(f'Rewrote {INDEX_MD}')


if __name__ == '__main__':
    main()
