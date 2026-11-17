# Servidor local mínimo para el Conversor 3D (solo escucha en este equipo).
# Sirve la carpeta "app" y abre el navegador. Cerrar esta ventana detiene el servidor.
$ErrorActionPreference = 'Stop'
$raiz = Join-Path $PSScriptRoot 'app'
$tipos = @{ '.html' = 'text/html; charset=utf-8'; '.js' = 'text/javascript; charset=utf-8'; '.mjs' = 'text/javascript; charset=utf-8';
  '.wasm' = 'application/wasm'; '.css' = 'text/css; charset=utf-8'; '.json' = 'application/json'; '.png' = 'image/png';
  '.svg' = 'image/svg+xml'; '.ico' = 'image/x-icon'; '.txt' = 'text/plain; charset=utf-8' }
$oyente = $null
foreach ($puerto in 8765..8790) {
  try {
    $oyente = New-Object System.Net.HttpListener
    $oyente.Prefixes.Add("http://localhost:$puerto/")
    $oyente.Start(); break
  } catch { $oyente = $null }
}
if (-not $oyente) { Write-Host 'No se pudo abrir ningún puerto local (8765-8790).'; Read-Host 'Pulsa Intro para salir'; exit 1 }
$url = "http://localhost:$puerto/"
Write-Host ''
Write-Host "  Conversor 3D en marcha: $url"
Write-Host '  Se abrirá en el navegador. No cierres esta ventana mientras lo uses.'
Write-Host '  Para terminar, cierra esta ventana.'
Write-Host ''
Start-Process $url
$raizCompleta = [System.IO.Path]::GetFullPath($raiz)
while ($oyente.IsListening) {
  $ctx = $oyente.GetContext()
  try {
    $ruta = [System.Uri]::UnescapeDataString($ctx.Request.Url.AbsolutePath.TrimStart('/'))
    if ($ruta -eq '') { $ruta = 'index.html' }
    $archivo = [System.IO.Path]::GetFullPath((Join-Path $raizCompleta $ruta))
    if (-not $archivo.StartsWith($raizCompleta) -or -not (Test-Path -LiteralPath $archivo -PathType Leaf)) {
      $ctx.Response.StatusCode = 404
    } else {
      $ext = [System.IO.Path]::GetExtension($archivo).ToLower()
      $ctx.Response.ContentType = if ($tipos.ContainsKey($ext)) { $tipos[$ext] } else { 'application/octet-stream' }
      $ctx.Response.Headers.Add('Cache-Control', 'no-cache')
      $datos = [System.IO.File]::ReadAllBytes($archivo)
      $ctx.Response.ContentLength64 = $datos.Length
      $ctx.Response.OutputStream.Write($datos, 0, $datos.Length)
    }
  } catch { try { $ctx.Response.StatusCode = 500 } catch {} }
  finally { try { $ctx.Response.OutputStream.Close() } catch {} }
}
