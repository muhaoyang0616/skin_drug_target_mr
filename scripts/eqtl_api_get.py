import sys
import time
import requests


def main():
    if len(sys.argv) != 4:
        raise SystemExit("usage: eqtl_api_get.py URL OUT TIMEOUT")
    url, out_path, timeout_s = sys.argv[1], sys.argv[2], float(sys.argv[3])
    if not url.startswith("https://www.ebi.ac.uk/eqtl/api/"):
        raise SystemExit("refusing non-eQTL-Catalogue URL")
    last_exc = None
    for attempt in range(4):
        try:
            response = requests.get(url, timeout=timeout_s)
            if 500 <= response.status_code < 600 and attempt < 3:
                time.sleep(2 * (attempt + 1))
                continue
            response.raise_for_status()
            with open(out_path, "wb") as f:
                f.write(response.content)
            return
        except requests.RequestException as exc:
            last_exc = exc
            if attempt < 3:
                time.sleep(2 * (attempt + 1))
                continue
            raise
    raise last_exc


if __name__ == "__main__":
    main()
