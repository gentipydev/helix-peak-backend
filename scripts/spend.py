"""What builds on demand have cost on Modal, and how near each cap is. Read-only.

    .venv/bin/python scripts/spend.py

From `spend_ledger` (0015), as the service and the resolver count it
(`app/budget.py`, `app/resolves.py`). Modal's own bill is `modal billing
summary`; this is the budget's view of it, which reserves a call's worst case
while it runs and settles what it cost when it ends.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import psycopg  # noqa: E402

from app import budget  # noqa: E402
from app.config import settings  # noqa: E402

_CALLS = """
select kind, count(*), coalesce(sum(usd), 0)::float8,
       count(*) filter (where finished_at is null)
from spend_ledger
where kind = any(%s::text[]) and at > now() - interval '1 day'
group by kind order by kind
"""
_ASKS = """
select count(*) filter (where kind = 'ask'), count(*) filter (where kind = 'resume'),
       count(distinct asker), max(n)
from spend_ledger
left join lateral (
    select count(*) as n from spend_ledger o
    where o.kind in ('ask', 'resume') and o.at > now() - interval '1 day'
      and o.asker is not distinct from spend_ledger.asker
) mine on true
where kind in ('ask', 'resume') and at > now() - interval '1 day'
"""
_HELD = """
select (select count(*) from protein where catalog_order is null),
       (select count(*) from resolve_request where state in ('queued', 'running'))
"""


def main() -> int:
    if not settings.database_url:
        raise SystemExit("DATABASE_URL is not set.")
    with psycopg.connect(settings.database_url, prepare_threshold=None) as conn:
        day, month = budget.spent(conn)
        settled_day, settled_month = budget.spent(conn, settled=True)
        calls = conn.execute(_CALLS, (list(budget.PAID),)).fetchall()
        asks, resumes, installs, most = conn.execute(_ASKS).fetchone()
        built, queued = conn.execute(_HELD).fetchone()
    print("Modal budget     last day ${:.4f} of ${:.2f}   this month ${:.4f} of ${:.2f}".format(
        day, budget.DAY_USD, month, budget.MONTH_USD))
    print("  settled only   last day ${:.4f}              this month ${:.4f}".format(
        settled_day, settled_month))
    for kind, count, usd, running in calls:
        print("  {:<9} {:>4} calls in the last day, ${:.4f}{}".format(
            kind, count, usd, ", {} running".format(running) if running else ""))
    print("Asks, last day   {} asks + {} resumes of {} (from {} installs; most from one: {} of {})"
          .format(asks, resumes, settings.resolves_per_day, installs, most or 0,
                  settings.resolves_per_install_per_day))
    print("Proteins held    {} built + {} queued of {}".format(
        built, queued, settings.max_built_proteins))
    print("A new build fits: {}".format(
        budget.short_of(settled_day, settled_month, budget.BUILD_USD) or "yes"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
