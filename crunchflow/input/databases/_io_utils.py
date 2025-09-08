from typing import List, Tuple

__all__ = [
    "_is_blank_or_comment",
    "_floats_in",
    "_tokens",
    "_separate_quoted_and_numeric",
    "_unquote",
    "_fmt_float",
    "_strip_bang",
]


# ----------------------------------------------------------------------------- #
# Utilities (local to this module)
# ----------------------------------------------------------------------------- #
def _strip_bang(s: str) -> str:
    """Remove trailing comments starting with '!'."""
    return s.split("!", 1)[0].strip()


def _is_blank_or_comment(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("!") or s.startswith("#")


def _floats_in(line: str) -> List[float]:
    vals: List[float] = []
    tok = ""
    for ch in line.replace(",", " ").strip():
        if ch in " \t":
            if tok:
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
                tok = ""
        else:
            tok += ch
    if tok:
        try:
            vals.append(float(tok))
        except ValueError:
            pass
    return vals


def _tokens(line: str) -> List[str]:
    # very light tokenizer: split on whitespace, keep quoted names intact
    out: List[str] = []
    buf: List[str] = []
    in_quote = False
    quote_char = ""
    for ch in line:
        if in_quote:
            buf.append(ch)
            if ch == quote_char:
                in_quote = False
                out.append("".join(buf).strip())
                buf = []
        else:
            if ch in ("'", '"'):
                in_quote = True
                quote_char = ch
                if buf:
                    if buf[-1].isspace():
                        buf = []
                    else:
                        out.append("".join(buf).strip())
                        buf = []
                buf = [ch]
            elif ch.isspace():
                if buf:
                    out.append("".join(buf).strip())
                    buf = []
            else:
                buf.append(ch)
    if buf:
        out.append("".join(buf).strip())
    return out


def _separate_quoted_and_numeric(toks: List[str]) -> Tuple[List[str], List[float]]:
    quoted: List[str] = []
    nums: List[float] = []
    for t in toks:
        if t.startswith("'") and t.endswith("'"):
            quoted.append(_unquote(t))
        else:
            try:
                nums.append(float(t.strip(",")))
            except ValueError:
                pass

    return quoted, nums


def _unquote(name: str) -> str:
    if (len(name) >= 2) and ((name[0] == "'" and name[-1] == "'") or (name[0] == '"' and name[-1] == '"')):
        return name[1:-1]
    return name


def _fmt_float(val: float, str_len: int = 10, align: str = "right") -> str:
    if (-1e-5 < val < 1e-5 and val != 0.0) or val > 1e8 or val < -1e8:
        # Split into base and exponent
        base, exp = f"{val:.8e}".split("e")
        base = base.rstrip("0")  # remove trailing zeros from base
        if base.endswith("."):
            base += "0"
        if align == "left":
            return f"{base}e{exp}".ljust(str_len)
        else:
            return f"{base}e{exp}".rjust(str_len)

    # Convert float to string, removing trailing zeros
    s = f"{val:.8f}".rstrip("0")
    # Add back in zero if we stripped all decimals
    if s.endswith("."):
        s += "0"

    if align == "left":
        return s.ljust(str_len)
    else:
        return s.rjust(str_len)
