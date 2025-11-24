#!/usr/bin/env bash
set -euo pipefail

# Robust downloader for ZINC curl command lists.
# Usage: ./zinc_retry_download.sh /path/to/ZINC22-downloader-3D-pdbqt.tgz.curl
# Env knobs:
#   ATTEMPTS (default 10)        : per-URL attempts before giving up
#   SLEEP_BASE (default 5)       : seconds between retries
#   BETWEEN (default 2)          : pause between different URLs
#   CONNECT_TIMEOUT (default 15) : curl connect timeout
#   MAX_TIME (default 900)       : max seconds per transfer

cmd_file="${1:-}"
if [[ -z "${cmd_file}" || ! -f "${cmd_file}" ]]; then
  echo "Usage: $0 path/to/commands.curl" >&2
  exit 1
fi

attempts="${ATTEMPTS:-10}"
sleep_base="${SLEEP_BASE:-5}"
between="${BETWEEN:-2}"
connect_timeout="${CONNECT_TIMEOUT:-15}"
max_time="${MAX_TIME:-900}"

success=0
fail=0
total=0

download_one() {
  local url="$1"
  local out="$2"
  mkdir -p "$(dirname "${out}")"

  local try
  for try in $(seq 1 "${attempts}"); do
    echo "[${try}/${attempts}] ${url} -> ${out}"
    if curl --fail --http1.1 --tlsv1.2 \
        --connect-timeout "${connect_timeout}" \
        --max-time "${max_time}" \
        --create-dirs -o "${out}.part" "${url}"; then
      mv -f "${out}.part" "${out}"
      echo "OK: ${out}"
      return 0
    fi
    rm -f "${out}.part"
    [[ "${try}" -lt "${attempts}" ]] && sleep "${sleep_base}"
  done
  return 1
}

while IFS= read -r line; do
  # Skip blanks/comments.
  [[ -z "${line// }" ]] && continue
  [[ "${line}" =~ ^# ]] && continue

  total=$((total + 1))
  out=""
  url=""
  prev=""

  for tok in ${line}; do
    if [[ "${prev}" == "-o" ]]; then
      out="${tok}"
    fi
    if [[ "${tok}" == http*://* ]]; then
      url="${tok}"
    fi
    prev="${tok}"
  done

  if [[ -z "${url}" || -z "${out}" ]]; then
    echo "Skip (unparsed): ${line}" >&2
    fail=$((fail + 1))
    continue
  fi

  if download_one "${url}" "${out}"; then
    success=$((success + 1))
  else
    echo "FAIL: ${url}" >&2
    fail=$((fail + 1))
  fi
  sleep "${between}"
done < "${cmd_file}"

echo "Done. success=${success} fail=${fail} total=${total}"
