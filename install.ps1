<#
.SYNOPSIS
    Install the isolf binary on Windows.

.DESCRIPTION
    Run it with:

        powershell -ExecutionPolicy ByPass -c "irm https://raw.githubusercontent.com/bondrewd/isolf/main/install.ps1 | iex"

    It downloads the prebuilt binary from the GitHub releases, checks its sha256,
    and installs it to %USERPROFILE%\.local\bin (added to your PATH).

    Pin a version or change the directory with the ISOLF_VERSION /
    ISOLF_INSTALL_DIR environment variables, or pass -Version / -InstallDir when
    running the downloaded file.
#>

param(
    [string]$Version = $env:ISOLF_VERSION,
    [string]$InstallDir = $env:ISOLF_INSTALL_DIR
)

$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'  # Invoke-WebRequest is much faster without it
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

$Repo = 'bondrewd/isolf'

function Say($m) { Write-Host "isolf: $m" }
function Fail($m) { Write-Host "isolf: error: $m" -ForegroundColor Red; exit 1 }

# Only x86_64 Windows binaries are published; they run on arm64 via emulation.
$target = 'x86_64-pc-windows-msvc'

# --- resolve the version -----------------------------------------------------
if ([string]::IsNullOrEmpty($Version)) {
    try {
        $latest = Invoke-RestMethod -Uri "https://api.github.com/repos/$Repo/releases/latest" -Headers @{ 'User-Agent' = 'isolf-install' }
        $tag = $latest.tag_name
    } catch {
        Fail "could not resolve the latest version (set ISOLF_VERSION to pin one)"
    }
} else {
    $tag = 'v' + ($Version -replace '^v', '')
}

$base = "isolf-$tag-$target"
$url = "https://github.com/$Repo/releases/download/$tag"

# --- download, verify, extract ----------------------------------------------
$tmp = Join-Path $env:TEMP ('isolf-' + [Guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $tmp | Out-Null
try {
    Say "downloading $base.zip"
    try {
        Invoke-WebRequest -UseBasicParsing -Uri "$url/$base.zip" -OutFile "$tmp\isolf.zip"
    } catch {
        Fail "download failed (does $tag exist?)"
    }

    try {
        Invoke-WebRequest -UseBasicParsing -Uri "$url/$base.sha256" -OutFile "$tmp\isolf.sha256"
        $expected = (((Get-Content "$tmp\isolf.sha256" -Raw) -split '\s+') | Where-Object { $_ })[0].ToLower()
        $got = (Get-FileHash "$tmp\isolf.zip" -Algorithm SHA256).Hash.ToLower()
        if ($expected -and ($expected -ne $got)) { Fail "checksum mismatch for $base.zip" }
    } catch {
        Say 'skipping checksum verification'
    }

    Expand-Archive -Path "$tmp\isolf.zip" -DestinationPath $tmp -Force
    $exe = Join-Path $tmp 'isolf.exe'
    if (-not (Test-Path $exe)) { Fail 'the archive did not contain isolf.exe' }

    # --- install -------------------------------------------------------------
    if ([string]::IsNullOrEmpty($InstallDir)) {
        $InstallDir = Join-Path $env:USERPROFILE '.local\bin'
    }
    New-Item -ItemType Directory -Force -Path $InstallDir | Out-Null
    Copy-Item -Force $exe (Join-Path $InstallDir 'isolf.exe')
    Say "installed isolf $($tag -replace '^v', '') to $InstallDir\isolf.exe"

    # --- add to the user PATH if missing -------------------------------------
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    if ($null -eq $userPath) { $userPath = '' }
    if (($userPath -split ';') -notcontains $InstallDir) {
        $newPath = if ($userPath) { "$userPath;$InstallDir" } else { $InstallDir }
        [Environment]::SetEnvironmentVariable('Path', $newPath, 'User')
        # Also update this session so isolf works right away, not only new terminals.
        if (($env:Path -split ';') -notcontains $InstallDir) { $env:Path = "$env:Path;$InstallDir" }
        Say "added $InstallDir to your PATH; run: isolf --version (new terminals pick it up automatically)"
    } else {
        Say 'run: isolf --version'
    }
} finally {
    Remove-Item -Recurse -Force $tmp -ErrorAction SilentlyContinue
}
