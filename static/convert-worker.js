// Conversión STL -> STEP en un Worker: reconoce superficies analíticas (planos,
// cilindros, conos, esferas y toros), escribe el STEP y lo valida con Open CASCADE.
import initOpenCascade from './vendor/replicad_single.js';
import { convertValidated } from './core/convert.js';

export const MAX_TRIANGLES = 400000;
const NAMES = { plane: ['plano', 'planos'], cylinder: ['cilindro', 'cilindros'], cone: ['cono', 'conos'], sphere: ['esfera', 'esferas'], torus: ['toro', 'toros'] };

function summary(report) {
  const parts = Object.keys(NAMES).filter(k => report.types[k]).map(k => `${report.types[k]} ${NAMES[k][report.types[k] === 1 ? 0 : 1]}`);
  let text = `${report.triangles.toLocaleString('es-ES')} triángulos → ${report.faces.toLocaleString('es-ES')} caras (${parts.join(', ')}).`;
  if (report.repairs) text += ' Se reparó la malla antes de convertir.';
  if (report.demoted) text += ' Algunas zonas se dejaron con caras planas.';
  return text;
}

self.onmessage = async ({ data }) => {
  try {
    const scale = { mm: 1, cm: 10, in: 25.4 }[data.units]; if (!scale) throw new Error('Unidad no válida.');
    const positions = new Float32Array(data.positions);
    if (!positions.length || positions.length % 9 || positions.length > MAX_TRIANGLES * 9) throw new Error(`STL vacío o superior a ${MAX_TRIANGLES.toLocaleString('es-ES')} triángulos.`);
    self.postMessage({ status: 'Cargando el motor CAD (22 MB la primera vez)…' });
    // El WebAssembly se compila una vez; cada motor nuevo reutiliza el módulo compilado.
    const wasm = await WebAssembly.compileStreaming(fetch(new URL('./vendor/replicad_single.wasm', import.meta.url)));
    const engine = () => initOpenCascade({ instantiateWasm: (imports, done) => { WebAssembly.instantiate(wasm, imports).then(i => done(i, wasm)); return {}; }, print: () => {}, printErr: () => {} });
    const result = await convertValidated(positions, engine, { scale, name: data.name || 'pieza', onStatus: status => self.postMessage({ status }) });
    const bytes = new TextEncoder().encode(result.step);
    self.postMessage({ result: bytes.buffer, summary: summary(result.report), report: result.report }, [bytes.buffer]);
  } catch (error) {
    self.postMessage({ error: typeof error === 'number' ? 'El motor CAD rechazó esta geometría.' : error.message || String(error) });
  }
};
