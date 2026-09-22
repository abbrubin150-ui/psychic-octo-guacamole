import * as SQLite from "expo-sqlite";
import type { Poem, PoemDraft, SortMode } from "./types";

type PoemRow = {
  id: number;
  title: string;
  body: string;
  written_at: string | null;
  created_at: string;
  updated_at: string;
  tags: string;
  favorite: number;
};

let dbPromise: Promise<SQLite.SQLiteDatabase> | null = null;
let migrated = false;

function rowToPoem(row: PoemRow): Poem {
  return {
    id: row.id,
    title: row.title,
    body: row.body,
    writtenAt: row.written_at,
    createdAt: row.created_at,
    updatedAt: row.updated_at,
    tags: row.tags,
    favorite: row.favorite === 1,
  };
}

async function migrate(db: SQLite.SQLiteDatabase) {
  if (migrated) return;

  await db.execAsync(`
    PRAGMA journal_mode = WAL;
    PRAGMA foreign_keys = ON;

    CREATE TABLE IF NOT EXISTS poems (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      title TEXT NOT NULL,
      body TEXT NOT NULL,
      written_at TEXT,
      created_at TEXT NOT NULL,
      updated_at TEXT NOT NULL,
      tags TEXT NOT NULL DEFAULT '',
      favorite INTEGER NOT NULL DEFAULT 0
    );

    CREATE INDEX IF NOT EXISTS idx_poems_written_at ON poems(written_at);
    CREATE INDEX IF NOT EXISTS idx_poems_updated_at ON poems(updated_at);
  `);

  migrated = true;
}

export async function getDb() {
  if (!dbPromise) dbPromise = SQLite.openDatabaseAsync("poemspace.db");
  const db = await dbPromise;
  await migrate(db);
  return db;
}

export async function listPoems(options?: {
  query?: string;
  month?: string;
  sort?: SortMode;
}): Promise<Poem[]> {
  const db = await getDb();
  const query = options?.query?.trim() ?? "";
  const month = options?.month?.trim() ?? "";
  const sort = options?.sort ?? "written_desc";

  const where: string[] = [];
  const params: (string | number | null)[] = [];

  if (query) {
    where.push("(title LIKE ? OR body LIKE ? OR tags LIKE ?)");
    const q = `%${query}%`;
    params.push(q, q, q);
  }

  if (/^\d{4}-\d{2}$/.test(month)) {
    where.push("written_at LIKE ?");
    params.push(`${month}%`);
  }

  const orderBy =
    sort === "written_asc"
      ? "COALESCE(written_at, '9999-12-31') ASC, id ASC"
      : sort === "updated_desc"
        ? "updated_at DESC, id DESC"
        : "COALESCE(written_at, '0000-00-00') DESC, id DESC";

  const sql = `
    SELECT id, title, body, written_at, created_at, updated_at, tags, favorite
    FROM poems
    ${where.length ? `WHERE ${where.join(" AND ")}` : ""}
    ORDER BY ${orderBy}
  `;

  const rows = await db.getAllAsync<PoemRow>(sql, params);
  return rows.map(rowToPoem);
}

export async function savePoem(draft: PoemDraft): Promise<number> {
  const db = await getDb();
  const now = new Date().toISOString();
  const title = draft.title.trim() || "ללא כותרת";
  const body = draft.body.trim();
  const tags = draft.tags.trim();
  const favorite = draft.favorite ? 1 : 0;

  if (draft.id) {
    await db.runAsync(
      `UPDATE poems
       SET title = ?, body = ?, written_at = ?, updated_at = ?, tags = ?, favorite = ?
       WHERE id = ?`,
      [title, body, draft.writtenAt, now, tags, favorite, draft.id],
    );
    return draft.id;
  }

  const result = await db.runAsync(
    `INSERT INTO poems
      (title, body, written_at, created_at, updated_at, tags, favorite)
     VALUES (?, ?, ?, ?, ?, ?, ?)`,
    [title, body, draft.writtenAt, now, now, tags, favorite],
  );

  return Number(result.lastInsertRowId);
}

export async function deletePoem(id: number) {
  const db = await getDb();
  await db.runAsync("DELETE FROM poems WHERE id = ?", [id]);
}

export async function toggleFavorite(id: number, favorite: boolean) {
  const db = await getDb();
  await db.runAsync(
    "UPDATE poems SET favorite = ?, updated_at = ? WHERE id = ?",
    [favorite ? 1 : 0, new Date().toISOString(), id],
  );
}
