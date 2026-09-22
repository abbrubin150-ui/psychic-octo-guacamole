import type { ClusterResult, Poem } from "./types";

const HEBREW_STOPWORDS = new Set([
  "אני", "את", "אתה", "אתם", "אתן", "הוא", "היא", "הם", "הן",
  "של", "שלי", "שלך", "שלו", "שלה", "שלנו", "שלכם", "שלהם",
  "עם", "על", "אל", "אם", "גם", "כי", "לא", "כן", "זה", "זו", "זאת",
  "מה", "מי", "כל", "עוד", "רק", "אבל", "או", "אז", "כמו", "היה", "הייתה",
  "יש", "אין", "בו", "בה", "לי", "לו", "לה", "לך", "לנו",
]);

function tokenize(text: string): string[] {
  return text
    .replace(/[\u0591-\u05C7]/g, "")
    .replace(/[^A-Za-z0-9\u0590-\u05FF\s]/g, " ")
    .toLowerCase()
    .split(/\s+/)
    .map((t) => t.trim())
    .filter((t) => t.length >= 2 && !HEBREW_STOPWORDS.has(t));
}

type SparseVector = Map<string, number>;

function buildTfIdf(poems: Poem[]) {
  const docs = poems.map((p) => tokenize(`${p.title} ${p.body}`));
  const df = new Map<string, number>();

  for (const tokens of docs) {
    for (const term of new Set(tokens)) {
      df.set(term, (df.get(term) ?? 0) + 1);
    }
  }

  const vectors: SparseVector[] = docs.map((tokens) => {
    const counts = new Map<string, number>();
    for (const term of tokens) counts.set(term, (counts.get(term) ?? 0) + 1);

    const vector = new Map<string, number>();
    const denom = Math.max(tokens.length, 1);

    for (const [term, count] of counts) {
      const idf = Math.log((poems.length + 1) / ((df.get(term) ?? 0) + 1)) + 1;
      vector.set(term, (count / denom) * idf);
    }

    let norm = 0;
    for (const value of vector.values()) norm += value * value;
    norm = Math.sqrt(norm) || 1;
    for (const [term, value] of vector) vector.set(term, value / norm);

    return vector;
  });

  return vectors;
}

function distance(a: SparseVector, b: SparseVector): number {
  let dot = 0;
  for (const [term, value] of a) dot += value * (b.get(term) ?? 0);
  return 1 - Math.max(-1, Math.min(1, dot));
}

function average(vectors: SparseVector[]): SparseVector {
  const out = new Map<string, number>();
  if (!vectors.length) return out;

  for (const vector of vectors) {
    for (const [term, value] of vector) {
      out.set(term, (out.get(term) ?? 0) + value / vectors.length);
    }
  }

  let norm = 0;
  for (const value of out.values()) norm += value * value;
  norm = Math.sqrt(norm) || 1;
  for (const [term, value] of out) out.set(term, value / norm);
  return out;
}

export function clusterPoems(poems: Poem[], requestedK: number): ClusterResult[] {
  if (!poems.length) return [];

  const vectors = buildTfIdf(poems);
  const k = Math.max(1, Math.min(Math.floor(requestedK), poems.length));

  const centroids: SparseVector[] = Array.from({ length: k }, (_, i) => {
    const index = Math.floor((i * poems.length) / k);
    return new Map(vectors[index] ?? []);
  });

  const assignments = new Array<number>(poems.length).fill(0);

  for (let iteration = 0; iteration < 30; iteration++) {
    let changed = false;

    for (let i = 0; i < vectors.length; i++) {
      const vector = vectors[i]!;
      let bestCluster = 0;
      let bestDistance = Number.POSITIVE_INFINITY;

      for (let c = 0; c < centroids.length; c++) {
        const d = distance(vector, centroids[c]!);
        if (d < bestDistance) {
          bestDistance = d;
          bestCluster = c;
        }
      }

      if (assignments[i] !== bestCluster) {
        assignments[i] = bestCluster;
        changed = true;
      }
    }

    const nextCentroids: SparseVector[] = [];
    for (let c = 0; c < k; c++) {
      const members = vectors.filter((_, i) => assignments[i] === c);
      nextCentroids.push(members.length ? average(members) : centroids[c]!);
    }

    centroids.splice(0, centroids.length, ...nextCentroids);
    if (!changed && iteration > 0) break;
  }

  return Array.from({ length: k }, (_, c) => {
    const poemIds = poems
      .filter((_, i) => assignments[i] === c)
      .map((poem) => poem.id);

    const topTerms = [...centroids[c]!.entries()]
      .sort((a, b) => b[1] - a[1])
      .slice(0, 8)
      .map(([term]) => term);

    return { id: c + 1, poemIds, topTerms };
  }).filter((cluster) => cluster.poemIds.length > 0);
}

export function corpusStats(poems: Poem[]) {
  const tokens = poems.flatMap((p) => tokenize(p.body));
  const unique = new Set(tokens);
  const dated = poems
    .map((p) => p.writtenAt)
    .filter((d): d is string => Boolean(d))
    .sort();

  const counts = new Map<string, number>();
  for (const token of tokens) counts.set(token, (counts.get(token) ?? 0) + 1);

  const topTerms = [...counts.entries()]
    .sort((a, b) => b[1] - a[1])
    .slice(0, 12)
    .map(([term, count]) => ({ term, count }));

  return {
    poems: poems.length,
    words: tokens.length,
    uniqueWords: unique.size,
    averageWords: poems.length ? Math.round(tokens.length / poems.length) : 0,
    firstDate: dated[0] ?? null,
    lastDate: dated[dated.length - 1] ?? null,
    topTerms,
  };
}
