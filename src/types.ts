export type Poem = {
  id: number;
  title: string;
  body: string;
  writtenAt: string | null;
  createdAt: string;
  updatedAt: string;
  tags: string;
  favorite: boolean;
};

export type PoemDraft = {
  id?: number;
  title: string;
  body: string;
  writtenAt: string | null;
  tags: string;
  favorite: boolean;
};

export type SortMode = "written_desc" | "written_asc" | "updated_desc";

export type ClusterResult = {
  id: number;
  poemIds: number[];
  topTerms: string[];
};
