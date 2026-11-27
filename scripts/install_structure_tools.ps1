$ErrorActionPreference = 'Stop'

# Downloads only the official command-line tools; binaries remain outside Git.
$Cache = Join-Path $env:USERPROFILE '.cache\geofuse-v2'
New-Item -ItemType Directory -Force -Path $Cache | Out-Null

$Assets = @(
    @{
        Name = 'AutoDock Vina 1.2.7'
        Url = 'https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe'
        Path = (Join-Path $Cache 'vina.exe')
        Sha256 = 'E0C4B2715E0C1A74F6E92D0F3BE0328AC97542EAFBC111E6B1EFAD897A73CCE5'
    },
    @{
        Name = 'P2Rank 2.5.1'
        Url = 'https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz'
        Path = (Join-Path $Cache 'p2rank_2.5.1.tar.gz')
        Sha256 = 'D243F2D9036AC053FEFB9407B5FE1C85F4FE077C519FD975AC585E995FEAB274'
    }
)

foreach ($Asset in $Assets) {
    if (!(Test-Path $Asset.Path)) {
        Write-Host "Downloading $($Asset.Name)..."
        Invoke-WebRequest -Uri $Asset.Url -OutFile $Asset.Path -TimeoutSec 600
    }
    $Actual = (Get-FileHash -Path $Asset.Path -Algorithm SHA256).Hash
    if ($Actual -ne $Asset.Sha256) {
        Remove-Item $Asset.Path -Force
        throw "SHA-256 verification failed for $($Asset.Name): $Actual"
    }
    Write-Host "Verified $($Asset.Name): $Actual"
}

$P2RankHome = Join-Path $Cache 'p2rank_2.5.1'
if (!(Test-Path (Join-Path $P2RankHome 'prank.bat'))) {
    tar -xzf (Join-Path $Cache 'p2rank_2.5.1.tar.gz') -C $Cache
}
if (!(Test-Path (Join-Path $P2RankHome 'prank.bat'))) {
    throw "P2Rank extraction failed: expected prank.bat under $P2RankHome"
}

$env:GEOFUSE_P2RANK_HOME = $P2RankHome
$env:GEOFUSE_VINA = Join-Path $Cache 'vina.exe'

Write-Host ''
Write-Host 'Structure tools are ready:'
Write-Host "  GEOFUSE_P2RANK_HOME=$env:GEOFUSE_P2RANK_HOME"
Write-Host "  GEOFUSE_VINA=$env:GEOFUSE_VINA"
Write-Host 'P2Rank requires Java 17 or later. Set JAVA_HOME to a JDK if it is not already configured.'
& $env:GEOFUSE_VINA --version
