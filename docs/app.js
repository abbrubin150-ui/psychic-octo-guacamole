// Utilities
function parsePvals(txt) {
  return txt.split(/[,\s]+/).map(s => parseFloat(s)).filter(x => !isNaN(x));
}

function sortWithIndex(arr) {
  return arr.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);
}

function unsort(valuesSorted, order) {
  const out = new Array(order.length);
  for (let k=0;k<order.length;k++) out[order[k]] = valuesSorted[k];
  return out;
}

// Holm step-down adjusted p-values
function holmAdjust(pvals) {
  const m = pvals.length;
  const ranked = sortWithIndex(pvals);
  const p = ranked.map(x=>x.v);
  let adj = new Array(m);
  let prev = 0;
  for (let j=0;j<m;j++){
    const factor = m - j;
    const val = Math.min(1, factor * p[j]);
    adj[j] = Math.max(prev, val);
    prev = adj[j];
  }
  const order = ranked.map(x=>x.i);
  return unsort(adj, order);
}

// BH adjusted p-values
function bhAdjust(pvals) {
  const m = pvals.length;
  const ranked = sortWithIndex(pvals);
  const p = ranked.map(x=>x.v);
  let adj = new Array(m);
  let next = 1;
  for (let j=m-1;j>=0;j--){
    const factor = m/(j+1);
    const val = Math.min(1, factor * p[j]);
    adj[j] = Math.min(next, val);
    next = adj[j];
  }
  const order = ranked.map(x=>x.i);
  return unsort(adj, order);
}

// BY adjusted p-values
function byAdjust(pvals) {
  const m = pvals.length;
  const c = Array.from({length:m}, (_,i)=>1/(i+1)).reduce((a,b)=>a+b,0);
  const ranked = sortWithIndex(pvals);
  const p = ranked.map(x=>x.v);
  let adj = new Array(m);
  let next = 1;
  for (let j=m-1;j>=0;j--){
    const factor = (m/(j+1)) * c;
    const val = Math.min(1, factor * p[j]);
    adj[j] = Math.min(next, val);
    next = adj[j];
  }
  const order = ranked.map(x=>x.i);
  return unsort(adj, order);
}

function run() {
  const ptxt = document.getElementById('pvals').value;
  const alpha = parseFloat(document.getElementById('alpha').value);
  const method = document.getElementById('method').value;
  const pvals = parsePvals(ptxt);
  if (!pvals.length) {
    document.getElementById('out').innerHTML = '<p>לא נמצאו ערכי p.</p>';
    return;
  }
  let padj;
  if (method === 'fdr_bh') padj = bhAdjust(pvals);
  else if (method === 'fdr_by') padj = byAdjust(pvals);
  else padj = holmAdjust(pvals);
  const reject = padj.map(v => v <= alpha);
  const rows = pvals.map((p,i)=>`<tr><td>${i+1}</td><td>${p.toFixed(6)}</td><td>${padj[i].toFixed(6)}</td><td>${reject[i] ? '<span class="badge">דחייה</span>' : ''}</td></tr>`).join('');
  document.getElementById('out').innerHTML = `
    <table>
      <thead><tr><th>#</th><th>p</th><th>p_adj</th><th>החלטה</th></tr></thead>
      <tbody>${rows}</tbody>
    </table>
    <p>שיטה: <b>${method}</b>, α/q=${alpha}</p>`;
}

document.getElementById('run').addEventListener('click', run);
