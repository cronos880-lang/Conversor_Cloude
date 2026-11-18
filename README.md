# Conversor 3D · STL → STEP con superficies reales

Versión derivada de [Conversor_3D](https://github.com/cronos880-lang/Conversor_3D). Usa la interfaz del Visor R21 y cambia por completo el motor de conversión: en lugar de crear una cara plana por cada triángulo, **reconoce las superficies de la pieza** (planos, cilindros, conos, esferas y toros) y escribe un STEP con esas caras, sus aristas (rectas, arcos y curvas) y sus vértices.

Todo sigue ocurriendo en el navegador: sin instalaciones, sin cuentas y sin subir archivos a ningún servidor.

## Versión descargable (sin conexión)

[Descargar Conversor_Cloude C5 (ZIP)](https://github.com/cronos880-lang/Conversor_Cloude/blob/claude/adoring-fermat-ptg1ro/entrega/Conversor_Cloude_C5.zip?raw=true)

Extrae todo el ZIP en una carpeta nueva y abre `Iniciar_Conversor.bat` (Windows). Arranca un servidor local mínimo con PowerShell, que ya viene con Windows, y abre el navegador; no instala nada ni necesita Python. En Linux o macOS usa `iniciar_conversor.sh` (requiere `python3`). Instrucciones en `LEEME_PRIMERO.txt` dentro del ZIP. El paquete se genera con `python tools/build_entrega.py C5`.

## Uso

1. Abre la web y selecciona tu STL.
2. Gira arrastrando, acerca con la rueda y desplaza con el botón derecho.
3. Elige las unidades y pulsa **Convertir y descargar STEP**.

Al terminar, la web indica cuántas caras tiene el STEP y de qué tipo son. El motor CAD (22 MB) se descarga la primera vez que se convierte.

## Qué hace la conversión

1. **Lee y repara la malla**: suelda vértices casi coincidentes, resuelve uniones en T, tapa huecos pequeños y admite piezas con cavidades internas o con varios cuerpos.
2. **Reconoce superficies**: agrupa los triángulos que pertenecen a un mismo plano, cilindro, cono, esfera o toro, con una tolerancia de 0,01 a 0,03 mm según el tamaño de la pieza para absorber el ruido habitual de los STL.
3. **Construye el sólido**: cada grupo pasa a ser una cara; los bordes entre caras se convierten en rectas, arcos de circunferencia o curvas spline, y los vértices se llevan sobre las superficies ajustadas.
4. **Valida el resultado** con Open CASCADE (WebAssembly): vuelve a abrir el STEP, comprueba que es un sólido cerrado y válido y que conserva volumen y área (margen del 2 %). Si alguna cara no supera la comprobación, se degrada solo esa zona: primero a contorno poligonal y, si hace falta, a caras planas. Nunca se descarga un STEP que no haya pasado la validación.

## Alcance real

- Funciona mejor con **piezas mecánicas** exportadas desde un programa CAD: taladros, ejes, chaflanes, redondeos sencillos.
- **No reconstruye el historial paramétrico ni las cotas.** El STEP es un sólido con caras reales, no un diseño con operaciones.
- Las **formas libres u orgánicas** (superficies esculpidas, perfiles spline, dentados, roscas) no se reconocen: esas zonas quedan con caras planas, igual que en una conversión facetada clásica.
- Los radios y ejes salen del ajuste a la malla; pueden diferir en centésimas del diseño original.
- Comprobado en Onshape (importación real por su API) con las piezas de prueba: maneta, Base y V_legoBrick4x2, que antes importaban con aviso, ya lo hacen sin él. No se garantiza para cualquier STL.
- Las relaciones de diseño evidentes (planos paralelos o perpendiculares, ejes coaxiales, redondeos tangentes) se imponen exactamente cuando la malla lo permite.

Límites: 100 MB y 400.000 triángulos. El STL no declara unidades: elige las correctas. El STEP se exporta en milímetros.

## Estructura

| Ruta | Contenido |
|---|---|
| `index.html`, `static/presentation.js` | Interfaz del Visor R21 (rama `perfeccion_Visualizador_R21` de Conversor_3D), sin cambios de diseño, conectada al motor de este repositorio. |
| `static/convert-worker.js` | Worker que ejecuta la conversión y la validación. |
| `static/core/mesh.js` | Lectura de STL, soldadura, reparación y adyacencias. |
| `static/core/fit.js` | Ajuste de cilindros, conos, esferas y toros. |
| `static/core/segment.js` | Reconocimiento de superficies sobre la malla. |
| `static/core/curves.js`, `brep.js` | Aristas, vértices y topología del sólido. |
| `static/core/step.js` | Escritura del STEP (AP214). |
| `static/core/validate.js`, `convert.js` | Validación con Open CASCADE y degradación de caras problemáticas. |
| `tests/` | Pruebas del núcleo y del flujo real en navegador. |
| `tools/` | Utilidades de desarrollo: conversión por lotes, diagnóstico e imágenes de la segmentación. |

## Desarrollo y pruebas

No hay compilación ni dependencias npm para ejecutar la aplicación. Sirve este directorio con cualquier servidor estático, por ejemplo `python -m http.server 8011`.

```bash
node tests/core_test.mjs
```

Prueba el núcleo con el motor WebAssembly real: caja, tubo, esfera, cavidad interna, varias piezas, unidades, reparación de huecos y rechazo de mallas no válidas.

```bash
node tests/browser_real.mjs carpeta_con_stl carpeta_salida
```

Convierte cada STL desde la página en un navegador (requiere Playwright y el servidor en el puerto 8011) y guarda los STEP descargados.

```bash
node tools/batch.mjs carpeta_con_stl carpeta_salida
```

Convierte por lotes en Node y vuelve a comprobar cada STEP con un Open CASCADE independiente (`pip install cadquery-ocp`).

## Publicación

Cada push a `main` (o a la rama de pruebas) ejecuta `.github/workflows/pages.yml`, que publica `index.html` y `static/` en GitHub Pages. Si GitHub no permite activarlo con el token de Actions, un administrador debe elegir **Settings → Pages → Source: GitHub Actions**.

## Componentes de terceros

- replicad-opencascadejs 1.1.0, licencia LGPL-2.1 incluida en `static/vendor/OPENCASCADE-LICENSE`. Distribución original: https://www.npmjs.com/package/replicad-opencascadejs . El JavaScript y el WebAssembly se conservan como archivos separados y reemplazables, sin modificaciones.
