"""What builds on demand may cost on Modal, and what they have cost.

One place, read by both sides: the service refuses a new ask while nothing is
left (`resolves.request`), and the resolver on Modal starts a call only where
its worst case still fits (`pipeline/resolver/modal_app.py`). Modal's images
carry `app/`, so both read these numbers; the service never imports `pipeline/`.

What is spent is `spend_ledger` (0015): a row for every call that costs money,
written before the call claims anything, with what the call could cost at most,
then replaced by what it did cost. A call that dies keeps its worst case. Its
sum over the last day and over this month (UTC, as Modal bills) is held to the
two limits below. Modal's own usage limit, set in its dashboard, stays the last
word: these keep the project well inside it, and say why a build waits.

What is not counted here is small and fixed: the sweep that runs every five
minutes, and the wake endpoint. Measured on 2026-11-13, a whole day of both was
about three cents.

Changing a limit is a code change: push it (the service) and deploy the
resolver again (`modal deploy pipeline/resolver/modal_app.py`).
"""

from typing import Optional, Tuple

# Dollars a rolling day of builds may cost on Modal, and dollars this
# calendar month (UTC). A typical build costs one to three cents (B2M, 102 s
# from ask to done on 2026-11-13); a scoring that runs to its timeout on the
# dearest GPU, about sixty-five. The month's is a third of the Starter plan's
# $30 of monthly credits, and half the $20 usage limit set on the workspace.
DAY_USD = 2.00
MONTH_USD = 10.00

# What one whole build may cost at most: a fold, its evidence and its scoring,
# each run to its timeout (`modal_app.WORST_USD`, which a test holds this to).
# The service takes an ask only while a build this dear still fits, so that a
# protein it accepts is one the resolver can afford to finish.
BUILD_USD = 0.85

# The calls on Modal that cost money, as `spend_ledger.kind` names them.
PAID = ("score", "fold", "evidence")

# The last day's and this month's spending. `settled` leaves out a call still
# running, which the service does not wait on; a call that never ended is
# counted at its worst case once an hour has passed, past every timeout.
_SPENT = """
select coalesce(sum(usd) filter (where at > now() - interval '1 day'), 0)::float8,
       coalesce(sum(usd) filter (where at >= date_trunc('month', now(), 'UTC')), 0)::float8
from spend_ledger
where kind = any(%s::text[])
  and at >= least(now() - interval '1 day', date_trunc('month', now(), 'UTC'))
  and (not %s or finished_at is not null or at < now() - interval '1 hour')
"""


def spent(conn, settled: bool = False) -> Tuple[float, float]:
    """Dollars the paid calls cost in the last day, and this month. A call
    still running counts at its worst case, unless `settled`."""
    row = conn.execute(_SPENT, (list(PAID), settled)).fetchone()
    return float(row[0]), float(row[1])


def short_of(day: float, month: float, need: float,
             day_usd: float = DAY_USD, month_usd: float = MONTH_USD) -> Optional[str]:
    """Why `need` more dollars do not fit, or None where they do."""
    if day + need > day_usd + 1e-9:
        return ("Builds on demand have used ${:.2f} of their ${:.2f} for the last day."
                .format(day, day_usd))
    if month + need > month_usd + 1e-9:
        return ("Builds on demand have used ${:.2f} of their ${:.2f} for this month."
                .format(month, month_usd))
    return None


# One reservation at a time, from reading what is spent to writing its row: the
# fold, the evidence and the scoring may each start at once.
_ONE_AT_A_TIME = "select pg_advisory_xact_lock(4848504542)"
_RESERVE = "insert into spend_ledger (kind, usd) values (%s, %s) returning id"
_SETTLE = """
update spend_ledger set usd = %s, seconds = %s, finished_at = now()
where id = %s and finished_at is null
"""


def reserve(conn, kind: str, usd: float) -> Tuple[Optional[int], Optional[str]]:
    """Set `usd`, a call's worst case, against the budget before it claims
    anything: its ledger row's id, or None and why it does not fit."""
    with conn.transaction():
        conn.execute(_ONE_AT_A_TIME)
        day, month = spent(conn)
        short = short_of(day, month, usd)
        if short is not None:
            return None, short
        return conn.execute(_RESERVE, (kind, round(usd, 6))).fetchone()[0], None


def settle(conn, run: int, usd: float, seconds: float) -> None:
    """Replace a call's worst case with what it cost, once it has ended."""
    conn.execute(_SETTLE, (round(usd, 6), seconds, run))
