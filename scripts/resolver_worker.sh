#!/bin/bash
# The resolver's worker on this Mac, as a launchd agent.
#
#   scripts/resolver_worker.sh install     write the agent and start it; it starts at each login
#   scripts/resolver_worker.sh uninstall   stop it and remove the agent; its log and state stay
#   scripts/resolver_worker.sh stop        stop it, and keep it stopped across logins
#   scripts/resolver_worker.sh start       let it run again after a stop
#   scripts/resolver_worker.sh restart     have it leave between two proteins and start again
#   scripts/resolver_worker.sh status      whether it runs, what waits in the queue, its last lines
#   scripts/resolver_worker.sh logs        follow its log
#
# The agent is pipeline/resolver/local_worker.py, kept running by launchd; the
# runbook is pipeline/resolver/README.md. The worker reads the repository's
# .env itself, so no value from it is written here, into the plist or anywhere.

set -euo pipefail

LABEL="com.helixpeek.resolver-worker"
BACKEND="$(cd "$(dirname "$0")/.." && pwd -P)"
PYTHON="$BACKEND/.worker-venv/bin/python"
ESM_PYTHON="$BACKEND/pipeline/.esm-venv/bin/python"
DOMAIN="gui/$(id -u)"
TARGET="$DOMAIN/$LABEL"
PLIST="$HOME/Library/LaunchAgents/$LABEL.plist"
LOGS="$HOME/Library/Logs/HelixPeek"
LOG="$LOGS/resolver-worker.log"
LAUNCHD_LOG="$LOGS/resolver-worker.launchd.log"
STATE="$HOME/Library/Application Support/HelixPeek/resolver-worker"
ME="scripts/resolver_worker.sh"

say() { printf '%s\n' "$*"; }
die() { printf '%s\n' "$*" >&2; exit 1; }

usage() { /usr/bin/sed -n '2,/^$/p' "$0" | /usr/bin/sed -e 's/^# \{0,1\}//'; }

loaded() { /bin/launchctl print "$TARGET" >/dev/null 2>&1; }

# The pid of the worker launchd is running, or nothing.
pid_of() {
    { /bin/launchctl print "$TARGET" 2>/dev/null || true; } \
        | /usr/bin/awk '/^\tpid = / && !found { print $3; found = 1 }'
}

disabled() {
    /bin/launchctl print-disabled "$DOMAIN" 2>/dev/null \
        | /usr/bin/grep -q -F "\"$LABEL\" => disabled"
}

# launchd sends SIGTERM and, if the worker has not left 30 s later
# (ExitTimeOut), SIGKILL. The worker leaves between two proteins, ending a
# scorer that is running so that its bake goes back on the queue.
unload() {
    loaded || return 0
    /bin/launchctl bootout "$TARGET" 2>/dev/null || true
    local tries=0
    while loaded; do
        tries=$((tries + 1))
        [ "$tries" -le 90 ] || die "launchd still holds $LABEL after 45 s."
        /bin/sleep 0.5
    done
}

require() {
    [ -x "$PYTHON" ] || die "There is no worker environment at .worker-venv. Make it with:
    /opt/homebrew/bin/python3.12 -m venv .worker-venv
    .worker-venv/bin/pip install -r requirements.txt"
    [ -x "$ESM_PYTHON" ] || die "There is no scorer environment at pipeline/.esm-venv.
pipeline/constraint/README.md says how it is made."
    [ -f "$BACKEND/.env" ] || die "There is no .env in $BACKEND. The worker reads
DATABASE_URL, SUPABASE_URL, SUPABASE_SERVICE_KEY and NCBI_EMAIL from it."
}

xml() { printf '%s' "$1" | /usr/bin/sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g'; }

# KeepAlive restarts the worker if it crashes, no sooner than ThrottleInterval
# after it last started, so a broken .env cannot spin. ProcessType is Standard
# because Background throttles CPU and I/O, which is the scoring. PATH is for
# caffeinate. stdout and stderr go to a file for what Python says before the
# worker's own log is open; the worker writes nothing else there.
write_plist() {
    /bin/mkdir -p "$(dirname "$PLIST")"
    local draft="$PLIST.$$"
    /bin/cat > "$draft" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>Label</key>
    <string>$LABEL</string>
    <key>ProgramArguments</key>
    <array>
        <string>$(xml "$PYTHON")</string>
        <string>-u</string>
        <string>-m</string>
        <string>pipeline.resolver.local_worker</string>
    </array>
    <key>WorkingDirectory</key>
    <string>$(xml "$BACKEND")</string>
    <key>RunAtLoad</key>
    <true/>
    <key>KeepAlive</key>
    <true/>
    <key>ThrottleInterval</key>
    <integer>60</integer>
    <key>ProcessType</key>
    <string>Standard</string>
    <key>ExitTimeOut</key>
    <integer>30</integer>
    <key>EnvironmentVariables</key>
    <dict>
        <key>PATH</key>
        <string>/usr/bin:/bin:/usr/sbin:/sbin</string>
        <key>PYTHONUNBUFFERED</key>
        <string>1</string>
    </dict>
    <key>StandardOutPath</key>
    <string>$(xml "$LAUNCHD_LOG")</string>
    <key>StandardErrorPath</key>
    <string>$(xml "$LAUNCHD_LOG")</string>
</dict>
</plist>
EOF
    if ! /usr/bin/plutil -lint "$draft" >/dev/null; then
        /bin/rm -f "$draft"
        die "The agent's plist did not pass plutil -lint; nothing was installed."
    fi
    /bin/chmod 644 "$draft"
    /bin/mv -f "$draft" "$PLIST"
}

do_install() {
    require
    /bin/mkdir -p "$LOGS" "$STATE"
    unload
    write_plist
    /bin/launchctl enable "$TARGET"
    /bin/launchctl bootstrap "$DOMAIN" "$PLIST"
    say "Installed the worker and started it. It starts again at each login."
    say "  agent  $PLIST"
    say "  log    $LOG"
    say "  state  $STATE"
}

do_uninstall() {
    unload
    if [ -f "$PLIST" ]; then
        /bin/rm -f "$PLIST"
        say "Stopped the worker and removed its agent."
    else
        say "No agent was installed."
    fi
    say "Its log stays in $LOGS, and its state in $STATE."
}

do_stop() {
    [ -f "$PLIST" ] || die "The agent is not installed, so there is nothing to stop."
    /bin/launchctl disable "$TARGET"
    if loaded; then
        unload
        say "Stopped the worker."
    else
        say "The worker was not running."
    fi
    say "It stays stopped across logins until: $ME start"
}

do_start() {
    [ -f "$PLIST" ] || die "The agent is not installed. Run: $ME install"
    require
    /bin/launchctl enable "$TARGET"
    if loaded; then
        /bin/launchctl kickstart "$TARGET" >/dev/null
        say "The worker was already loaded, and is running (pid $(pid_of))."
    else
        /bin/launchctl bootstrap "$DOMAIN" "$PLIST"
        say "Started the worker. It starts again at each login."
    fi
}

# SIGTERM, which the worker answers by leaving between two proteins, and then a
# kickstart. launchd's own stop (`kickstart -k`, `bootout`) sends SIGTERM too,
# and SIGKILL once ExitTimeOut is up, which a slow download would not survive.
# launchd starts no process sooner than ThrottleInterval after the last one,
# and the kickstart waits for that: within a minute of the worker's last start,
# a restart takes up to a minute to return.
do_restart() {
    loaded || die "The worker is not loaded. Run: $ME start"
    local old tries=0
    old="$(pid_of)"
    if [ -n "$old" ]; then
        /bin/launchctl kill SIGTERM "$TARGET"
        while [ "$(pid_of)" = "$old" ]; do
            tries=$((tries + 1))
            if [ "$tries" -gt 120 ]; then
                say "The worker (pid $old) has been asked to leave and is finishing the protein in hand."
                say "launchd starts it again when it has."
                return 0
            fi
            /bin/sleep 0.5
        done
    fi
    /bin/launchctl kickstart "$TARGET" >/dev/null
    say "Restarted the worker (pid $(pid_of))."
}

do_status() {
    if loaded; then
        /bin/launchctl print "$TARGET" \
            | /usr/bin/grep -E $'^\t(state|pid|runs|last exit code) = ' \
            | /usr/bin/sed -e $'s/^\t//' || true
    elif [ ! -f "$PLIST" ]; then
        say "not installed ($ME install)"
    elif disabled; then
        say "stopped, and staying stopped across logins ($ME start)"
    else
        say "installed, not loaded: it loads at the next login ($ME start)"
    fi
    say ""
    if [ -x "$PYTHON" ]; then
        (cd "$BACKEND" && "$PYTHON" -u -m pipeline.resolver.local_worker --queue) || true
    else
        say "the queue was not read: there is no .worker-venv"
    fi
    say ""
    if [ -f "$LOG" ]; then
        say "$LOG"
        /usr/bin/tail -n 10 "$LOG" | /usr/bin/sed -e 's/^/  /'
    else
        say "no log yet at $LOG"
    fi
    if [ -s "$LAUNCHD_LOG" ]; then
        say ""
        say "$LAUNCHD_LOG (what was said outside the log)"
        /usr/bin/tail -n 5 "$LAUNCHD_LOG" | /usr/bin/sed -e 's/^/  /'
    fi
}

do_logs() {
    [ -f "$LOG" ] || die "There is no log yet at $LOG."
    # -F: the log is rotated at 1 MB, and the file to follow is then a new one.
    exec /usr/bin/tail -n 50 -F "$LOG"
}

case "${1:-}" in
    install) do_install ;;
    uninstall) do_uninstall ;;
    stop) do_stop ;;
    start) do_start ;;
    restart) do_restart ;;
    status) do_status ;;
    logs) do_logs ;;
    *) usage; exit 2 ;;
esac
