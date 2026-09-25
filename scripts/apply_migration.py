"""Apply one migration file to the database in `.env`, as `psql -f` would.

    .venv/bin/python scripts/apply_migration.py migrations/0002_protein_index.sql

For a machine without psql. Each migration carries its own `begin`/`commit`,
so the file is sent as it is on an autocommit connection: psycopg must not
wrap it in a transaction of its own, or the file's `commit` would end that one
early.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import psycopg  # noqa: E402

from app.config import settings  # noqa: E402


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    if not settings.database_url:
        raise SystemExit("DATABASE_URL is not set.")
    with open(sys.argv[1]) as handle:
        sql = handle.read()
    with psycopg.connect(settings.database_url, autocommit=True, prepare_threshold=None) as conn:
        conn.execute(sql)
    print("Applied {}.".format(sys.argv[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
