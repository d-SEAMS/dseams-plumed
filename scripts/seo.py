#!/usr/bin/env python3
"""Add and validate SEO metadata for every generated HTML page."""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from pathlib import Path
from urllib.parse import urljoin
from xml.etree import ElementTree

SITE_NAME = "dseams-plumed"
SITE_DESCRIPTION = (
    "PLUMED action for d-SEAMS topological ice and cage analysis with the Chemical "
    "Identification of Local Liquid-like and Solid-like Habitats Plus classifier"
)
SITE_URL = "https://plumed.dseams.info/"
REPOSITORY = "https://github.com/d-SEAMS/dseams-plumed"
IMAGE = "https://dseams.info/images/og-card.png"
DOI = "https://doi.org/10.1021/acs.jcim.0c00031"

TITLE_RE = re.compile(r"<title\b[^>]*>(.*?)</title>", re.I | re.S)
HEAD_END_RE = re.compile(r"</head\s*>", re.I)
SEO_RE = [
    re.compile(r'<meta\b[^>]*(?:name|property)=["\'](?:description|og:[^"\']+|twitter:[^"\']+)["\'][^>]*>\s*', re.I),
    re.compile(r'<link\b[^>]*rel=["\']canonical["\'][^>]*>\s*', re.I),
    re.compile(r'<script\b[^>]*id=["\']seo-json-ld["\'][^>]*>.*?</script>\s*', re.I | re.S),
]


def pages(root: Path) -> list[Path]:
    return sorted(root.rglob("*.html"))


def canonical_url(root: Path, page: Path) -> str:
    rel = page.relative_to(root)
    if rel.name == "index.html":
        path = "" if rel.parent == Path(".") else rel.parent.as_posix() + "/"
    else:
        path = rel.as_posix()
    return urljoin(SITE_URL, path)


def page_title(text: str, page: Path, root: Path) -> str:
    rel = page.relative_to(root).as_posix()
    if rel == "404.html":
        return f"Page not found - {SITE_NAME}"
    match = TITLE_RE.search(text)
    title = re.sub(r"<[^>]+>", "", match.group(1)) if match else ""
    title = " ".join(html.unescape(title).split())
    if not title:
        title = page.relative_to(root).with_suffix("").as_posix().replace("/", " - ")
    return title


def is_error_page(root: Path, page: Path) -> bool:
    return page.relative_to(root).as_posix() == "404.html"


def metadata(title: str, description: str, canonical: str) -> str:
    attr = lambda value: html.escape(value, quote=True)
    structured = {
        "@context": "https://schema.org",
        "@graph": [
            {
                "@type": "SoftwareApplication",
                "name": SITE_NAME,
                "description": SITE_DESCRIPTION,
                "url": SITE_URL,
                "applicationCategory": "DeveloperApplication",
                "citation": DOI,
            },
            {
                "@type": "SoftwareSourceCode",
                "name": SITE_NAME,
                "description": SITE_DESCRIPTION,
                "url": SITE_URL,
                "codeRepository": REPOSITORY,
                "programmingLanguage": ["C++", "PLUMED"],
                "license": "https://opensource.org/license/mit",
                "citation": DOI,
            },
        ],
    }
    tags = [
        f'<meta name="description" content="{attr(description)}">',
        f'<link rel="canonical" href="{attr(canonical)}">',
        '<meta property="og:type" content="website">',
        f'<meta property="og:site_name" content="{attr(SITE_NAME)}">',
        f'<meta property="og:title" content="{attr(title)}">',
        f'<meta property="og:description" content="{attr(description)}">',
        f'<meta property="og:url" content="{attr(canonical)}">',
        f'<meta property="og:image" content="{attr(IMAGE)}">',
        '<meta name="twitter:card" content="summary_large_image">',
        f'<meta name="twitter:title" content="{attr(title)}">',
        f'<meta name="twitter:description" content="{attr(description)}">',
        f'<meta name="twitter:image" content="{attr(IMAGE)}">',
        '<script id="seo-json-ld" type="application/ld+json">'
        + json.dumps(structured, separators=(",", ":"), ensure_ascii=True)
        + "</script>",
    ]
    return "\n".join(tags)


def enhance_site(root: Path) -> None:
    root = root.resolve()
    html_pages = pages(root)
    if not html_pages:
        raise RuntimeError(f"no HTML pages found in {root}")
    seen: dict[str, Path] = {}
    for page in html_pages:
        text = page.read_text(encoding="utf-8")
        title = page_title(text, page, root)
        if title in seen:
            label = page.relative_to(root).with_suffix("").as_posix().replace("/", " - ")
            title = f"{title} ({label})"
        seen[title] = page
        description = f"{title}. {SITE_DESCRIPTION}."
        canonical = canonical_url(root, page)
        text = TITLE_RE.sub(f"<title>{html.escape(title)}</title>", text, count=1)
        for pattern in SEO_RE:
            text = pattern.sub("", text)
        if not HEAD_END_RE.search(text):
            raise RuntimeError(f"{page} has no closing head element")
        text = HEAD_END_RE.sub(metadata(title, description, canonical) + "\n</head>", text, count=1)
        page.write_text(text, encoding="utf-8")
    sitemap = ElementTree.Element(
        "urlset", xmlns="http://www.sitemaps.org/schemas/sitemap/0.9"
    )
    for page in html_pages:
        if is_error_page(root, page):
            continue
        entry = ElementTree.SubElement(sitemap, "url")
        ElementTree.SubElement(entry, "loc").text = canonical_url(root, page)
    ElementTree.ElementTree(sitemap).write(
        root / "sitemap.xml", encoding="utf-8", xml_declaration=True
    )
    (root / "robots.txt").write_text(
        f"User-agent: *\nAllow: /\nSitemap: {urljoin(SITE_URL, 'sitemap.xml')}\n",
        encoding="utf-8",
    )


def check_site(root: Path) -> None:
    root = root.resolve()
    html_pages = pages(root)
    if not html_pages:
        raise RuntimeError(f"no HTML pages found in {root}")
    titles: set[str] = set()
    descriptions: set[str] = set()
    canonicals: set[str] = set()
    sitemap_canonicals: set[str] = set()
    required = {
        ("name", "description"),
        ("property", "og:type"),
        ("property", "og:site_name"),
        ("property", "og:title"),
        ("property", "og:description"),
        ("property", "og:url"),
        ("property", "og:image"),
        ("name", "twitter:card"),
        ("name", "twitter:title"),
        ("name", "twitter:description"),
        ("name", "twitter:image"),
    }
    for page in html_pages:
        text = page.read_text(encoding="utf-8")
        title_matches = TITLE_RE.findall(text)
        if len(title_matches) != 1:
            raise RuntimeError(f"{page}: expected exactly one title")
        title = " ".join(html.unescape(title_matches[0]).split())
        tags = {
            (kind.lower(), name.lower()): value
            for kind, name, value in re.findall(
                r'<meta\b[^>]*(name|property)=["\']([^"\']+)["\'][^>]*content=["\']([^"\']*)["\'][^>]*>',
                text,
                re.I,
            )
        }
        missing = required - tags.keys()
        if missing:
            raise RuntimeError(f"{page}: missing metadata {sorted(missing)}")
        og_images = re.findall(
            r'<meta\b[^>]*property=["\']og:image["\'][^>]*content=["\']([^"\']*)["\'][^>]*>',
            text,
            re.I,
        )
        if og_images != [IMAGE]:
            raise RuntimeError(f"{page}: expected exactly one non-empty og:image")
        if tags[("name", "twitter:image")] != IMAGE:
            raise RuntimeError(f"{page}: twitter:image must be {IMAGE}")
        description = html.unescape(tags[("name", "description")])
        canonical_matches = re.findall(
            r'<link\b[^>]*rel=["\']canonical["\'][^>]*href=["\']([^"\']+)["\'][^>]*>',
            text,
            re.I,
        )
        expected = canonical_url(root, page)
        if canonical_matches != [expected]:
            raise RuntimeError(f"{page}: canonical is not {expected}")
        if tags[("property", "og:url")] != expected:
            raise RuntimeError(f"{page}: og:url does not match its canonical")
        if title in titles or description in descriptions or expected in canonicals:
            raise RuntimeError(f"{page}: title, description, and canonical must be unique")
        if page.relative_to(root).as_posix() == "index.html" and title != SITE_NAME:
            raise RuntimeError(f"{page}: home title must be {SITE_NAME!r}, got {title!r}")
        titles.add(title)
        descriptions.add(description)
        canonicals.add(expected)
        if not is_error_page(root, page):
            sitemap_canonicals.add(expected)
        scripts = re.findall(
            r'<script\b[^>]*type=["\']application/ld\+json["\'][^>]*>(.*?)</script>',
            text,
            re.I | re.S,
        )
        if not any(
            {item.get("@type") for item in json.loads(script).get("@graph", [])}
            >= {"SoftwareApplication", "SoftwareSourceCode"}
            and DOI in script
            for script in scripts
        ):
            raise RuntimeError(f"{page}: missing required software JSON-LD")
    sitemap = ElementTree.parse(root / "sitemap.xml")
    locations = {node.text for node in sitemap.findall(".//{*}loc")}
    if locations != sitemap_canonicals:
        raise RuntimeError("sitemap URLs do not match generated HTML canonicals")
    robots = (root / "robots.txt").read_text(encoding="utf-8").splitlines()
    if robots != [
        "User-agent: *",
        "Allow: /",
        f"Sitemap: {urljoin(SITE_URL, 'sitemap.xml')}",
    ]:
        raise RuntimeError("robots.txt must allow crawling and name the custom-domain sitemap")
    print(f"SEO validated for {len(html_pages)} HTML pages")


def on_post_build(config, **kwargs) -> None:
    enhance_site(Path(config["site_dir"]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("site", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    (check_site if args.check else enhance_site)(args.site)
    return 0


if __name__ == "__main__":
    sys.exit(main())
