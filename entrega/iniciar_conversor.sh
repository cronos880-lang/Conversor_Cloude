#!/bin/sh
# Linux / macOS: sirve la carpeta "app" en este equipo y abre el navegador.
cd "$(dirname "$0")/app" || exit 1
PUERTO=8765
URL="http://localhost:$PUERTO/"
( sleep 1; (xdg-open "$URL" || open "$URL") >/dev/null 2>&1 ) &
echo "Conversor 3D en marcha: $URL (Ctrl+C para terminar)"
exec python3 -m http.server "$PUERTO" --bind 127.0.0.1
