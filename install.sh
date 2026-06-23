#!/bin/sh
# Install the isolf binary on macOS or Linux.
#
#   curl -LsSf https://raw.githubusercontent.com/bondrewd/isolf/main/install.sh | sh
#   wget -qO-  https://raw.githubusercontent.com/bondrewd/isolf/main/install.sh | sh
#
# It downloads the right prebuilt binary from the GitHub releases, checks its
# sha256, and installs it to ~/.local/bin (a `isolf --version` away from done).
#
# Knobs (environment variables):
#   ISOLF_VERSION       install a specific version, e.g. 0.1.1 (default: latest)
#   ISOLF_INSTALL_DIR   install directory (default: $XDG_BIN_HOME or ~/.local/bin)
# A version may also be passed as the first argument: `... | sh -s -- 0.1.1`.

set -eu

REPO="bondrewd/isolf"

say() { printf 'isolf: %s\n' "$1"; }
err() { printf 'isolf: error: %s\n' "$1" >&2; exit 1; }

# --- pick a downloader -------------------------------------------------------
if command -v curl >/dev/null 2>&1; then
    DL="curl"
elif command -v wget >/dev/null 2>&1; then
    DL="wget"
else
    err "this installer needs curl or wget"
fi

# fetch <url>            -> body on stdout
fetch() {
    if [ "$DL" = curl ]; then curl -fsSL "$1"; else wget -qO- "$1"; fi
}
# download <url> <dest>  -> save to file
download() {
    if [ "$DL" = curl ]; then curl -fsSL -o "$2" "$1"; else wget -qO "$2" "$1"; fi
}

# --- detect the platform -----------------------------------------------------
os="$(uname -s)"
arch="$(uname -m)"
case "$os" in
    Darwin) target="universal-apple-darwin" ;;
    Linux)
        case "$arch" in
            x86_64 | amd64) target="x86_64-unknown-linux-musl" ;;
            aarch64 | arm64) target="aarch64-unknown-linux-musl" ;;
            *) err "unsupported Linux architecture '$arch'; see https://github.com/$REPO/releases" ;;
        esac
        ;;
    *) err "unsupported OS '$os'; on Windows download the .zip from https://github.com/$REPO/releases" ;;
esac

# --- resolve the version -----------------------------------------------------
version="${1:-${ISOLF_VERSION:-}}"
if [ -n "$version" ]; then
    tag="v${version#v}"
else
    tag="$(fetch "https://api.github.com/repos/$REPO/releases/latest" \
        | grep -m1 '"tag_name"' \
        | sed -E 's/.*"tag_name" *: *"([^"]+)".*/\1/')"
    [ -n "$tag" ] || err "could not resolve the latest version (set ISOLF_VERSION to pin one)"
fi

base="isolf-$tag-$target"
url="https://github.com/$REPO/releases/download/$tag"

# --- download, verify, extract ----------------------------------------------
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT INT TERM

say "downloading $base.tar.gz"
download "$url/$base.tar.gz" "$tmp/archive.tar.gz" || err "download failed (does $tag exist?)"

if download "$url/$base.sha256" "$tmp/archive.sha256" 2>/dev/null; then
    expected="$(awk '{print $1; exit}' "$tmp/archive.sha256")"
    if command -v sha256sum >/dev/null 2>&1; then
        got="$(sha256sum "$tmp/archive.tar.gz" | awk '{print $1}')"
    elif command -v shasum >/dev/null 2>&1; then
        got="$(shasum -a 256 "$tmp/archive.tar.gz" | awk '{print $1}')"
    else
        got=""
        say "no sha256 tool found; skipping checksum verification"
    fi
    [ -z "$got" ] || [ "$got" = "$expected" ] || err "checksum mismatch for $base.tar.gz"
fi

tar -xzf "$tmp/archive.tar.gz" -C "$tmp"
[ -f "$tmp/isolf" ] || err "the archive did not contain an isolf binary"

# --- install -----------------------------------------------------------------
install_dir="${ISOLF_INSTALL_DIR:-${XDG_BIN_HOME:-$HOME/.local/bin}}"
mkdir -p "$install_dir"
if ! install -m 0755 "$tmp/isolf" "$install_dir/isolf" 2>/dev/null; then
    cp "$tmp/isolf" "$install_dir/isolf"
    chmod 0755 "$install_dir/isolf"
fi
# A binary fetched by a script carries no quarantine flag, but clear it anyway in
# case the directory was flagged, so the first run is not blocked on macOS.
[ "$os" = Darwin ] && xattr -d com.apple.quarantine "$install_dir/isolf" 2>/dev/null || true

say "installed isolf ${tag#v} to $install_dir/isolf"

# --- PATH guidance -----------------------------------------------------------
# If the install directory is already on PATH, isolf is ready. Otherwise show how
# to add it, tailored to the user's login shell — the installer runs under
# /bin/sh, so $SHELL (the login shell) is the signal to use, not $0.
case ":$PATH:" in
    *":$install_dir:"*)
        say "run: isolf --version"
        ;;
    *)
        # Keep $PATH literal so the printed line is copy-paste correct.
        export_line="export PATH=\"$install_dir:\$PATH\""
        say "$install_dir is not on your PATH."
        case "${SHELL##*/}" in
            fish)
                say "add it (fish):"
                printf '\n    fish_add_path %s\n\n' "$install_dir"
                ;;
            zsh)
                say "add it to ~/.zshrc (zsh):"
                printf "\n    echo '%s' >> ~/.zshrc && source ~/.zshrc\n\n" "$export_line"
                ;;
            bash)
                say "add it to ~/.bashrc (bash):"
                printf "\n    echo '%s' >> ~/.bashrc && source ~/.bashrc\n\n" "$export_line"
                ;;
            *)
                say "add it to your shell startup file (e.g. ~/.profile):"
                printf '\n    %s\n\n' "$export_line"
                ;;
        esac
        say "then run: isolf --version"
        ;;
esac
