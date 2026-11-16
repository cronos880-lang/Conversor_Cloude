import fs from 'fs';
const p = JSON.parse(fs.readFileSync(process.argv[2])); const self = new Map(); const dt = p.timeDeltas; const idx = new Map(p.nodes.map(n => [n.id, n]));
p.samples.forEach((id, i) => { const n = idx.get(id); const k = `${n.callFrame.functionName || '(anon)'}:${n.callFrame.lineNumber + 1} ${n.callFrame.url.split('/').pop()}`; self.set(k, (self.get(k) || 0) + dt[i]); });
const tot = [...self.values()].reduce((a, b) => a + b, 0);
[...self.entries()].sort((a, b) => b[1] - a[1]).slice(0, 18).forEach(([k, v]) => console.log((v / tot * 100).toFixed(1).padStart(5) + '% ' + k));
