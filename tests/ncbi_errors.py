"""Real NCBI failure bodies, captured so the tests exercise the actual shapes.

NCBI does not use status codes consistently and does not signal "no such
record" with an empty body. The spaced-out letters are literal: NCBI really
sends "F a i l e d".
"""

import email.message
import io
from urllib.error import HTTPError

UNPARSEABLE_ID_BODY = (
    "Error: F a i l e d  t o  u n d e r s t a n d  i d :  N O T _ A _ R E A L _ I D \n\n\n"
)
MISSING_ACCESSION_BODY = (
    b"Error: CEFetchPApplication::proxy_stream(): Error: "
    b"F a i l e d  t o  r e t r i e v e  s e q u e n c e :  N G _ 9 9 9 9 9 9 \n\n"
)


def http_error(code, body):
    """Build a real HTTPError, so exc.read() works as it does in production."""
    return HTTPError(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        code=code,
        msg="Bad Request",
        hdrs=email.message.Message(),
        fp=io.BytesIO(body),
    )
