import { StatusBar } from "expo-status-bar";
import React, { useCallback, useEffect, useMemo, useState } from "react";
import {
  Alert,
  FlatList,
  I18nManager,
  KeyboardAvoidingView,
  Platform,
  Pressable,
  SafeAreaView,
  ScrollView,
  StyleSheet,
  Text,
  TextInput,
  View,
} from "react-native";
import { clusterPoems, corpusStats } from "./src/analysis";
import { deletePoem, listPoems, savePoem, toggleFavorite } from "./src/db";
import { win95 } from "./src/theme";
import type { ClusterResult, Poem, PoemDraft, SortMode } from "./src/types";

I18nManager.allowRTL(true);

type Screen = "archive" | "editor" | "clusters" | "stats";

const emptyDraft = (): PoemDraft => ({
  title: "",
  body: "",
  writtenAt: new Date().toISOString().slice(0, 10),
  tags: "",
  favorite: false,
});

function RetroButton({
  children,
  onPress,
  active = false,
  disabled = false,
}: {
  children: React.ReactNode;
  onPress: () => void;
  active?: boolean;
  disabled?: boolean;
}) {
  return (
    <Pressable
      accessibilityRole="button"
      disabled={disabled}
      onPress={onPress}
      style={({ pressed }) => [
        styles.button,
        active && styles.buttonActive,
        pressed && !active && styles.buttonPressed,
        disabled && styles.disabled,
      ]}
    >
      <Text style={[styles.buttonText, active && styles.buttonTextActive]}>
        {children}
      </Text>
    </Pressable>
  );
}

function RetroField({
  value,
  onChangeText,
  placeholder,
  multiline = false,
}: {
  value: string;
  onChangeText: (value: string) => void;
  placeholder?: string;
  multiline?: boolean;
}) {
  return (
    <TextInput
      value={value}
      onChangeText={onChangeText}
      placeholder={placeholder}
      placeholderTextColor={win95.mid}
      multiline={multiline}
      textAlign="right"
      style={[styles.field, multiline && styles.textArea]}
    />
  );
}

function Window({
  title,
  children,
}: {
  title: string;
  children: React.ReactNode;
}) {
  return (
    <View style={styles.window}>
      <View style={styles.titleBar}>
        <Text style={styles.titleBarText}>{title}</Text>
        <View style={styles.titleButtons}>
          <Text style={styles.titleMini}>_</Text>
          <Text style={styles.titleMini}>□</Text>
          <Text style={styles.titleMini}>×</Text>
        </View>
      </View>
      <View style={styles.windowBody}>{children}</View>
    </View>
  );
}

export default function App() {
  const [screen, setScreen] = useState<Screen>("archive");
  const [poems, setPoems] = useState<Poem[]>([]);
  const [corpus, setCorpus] = useState<Poem[]>([]);
  const [query, setQuery] = useState("");
  const [month, setMonth] = useState("");
  const [sort, setSort] = useState<SortMode>("written_desc");
  const [draft, setDraft] = useState<PoemDraft>(emptyDraft());
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
  const [k, setK] = useState(3);
  const [clusters, setClusters] = useState<ClusterResult[]>([]);

  const refresh = useCallback(async () => {
    try {
      setError("");
      const [filtered, all] = await Promise.all([
        listPoems({ query, month, sort }),
        listPoems({ sort: "written_desc" }),
      ]);
      setPoems(filtered);
      setCorpus(all);
    } catch (e) {
      setError(e instanceof Error ? e.message : "שגיאה לא ידועה");
    }
  }, [query, month, sort]);

  useEffect(() => {
    void refresh();
  }, [refresh]);

  const stats = useMemo(() => corpusStats(corpus), [corpus]);

  function editPoem(poem: Poem) {
    setDraft({
      id: poem.id,
      title: poem.title,
      body: poem.body,
      writtenAt: poem.writtenAt,
      tags: poem.tags,
      favorite: poem.favorite,
    });
    setScreen("editor");
  }

  function newPoem() {
    setDraft(emptyDraft());
    setScreen("editor");
  }

  async function persistDraft() {
    if (!draft.body.trim()) {
      Alert.alert("PoemSpace 95", "צריך לכתוב טקסט לפני השמירה.");
      return;
    }
    if (draft.writtenAt && !/^\d{4}-\d{2}-\d{2}$/.test(draft.writtenAt)) {
      Alert.alert("PoemSpace 95", "תאריך צריך להיות בפורמט YYYY-MM-DD.");
      return;
    }

    try {
      setBusy(true);
      await savePoem(draft);
      setDraft(emptyDraft());
      await refresh();
      setScreen("archive");
    } catch (e) {
      Alert.alert("שגיאת שמירה", e instanceof Error ? e.message : "לא ניתן לשמור");
    } finally {
      setBusy(false);
    }
  }

  function confirmDelete(poem: Poem) {
    Alert.alert(
      "מחיקת שיר",
      `למחוק לצמיתות את "${poem.title}"?`,
      [
        { text: "ביטול", style: "cancel" },
        {
          text: "מחק",
          style: "destructive",
          onPress: () => {
            void (async () => {
              await deletePoem(poem.id);
              await refresh();
            })();
          },
        },
      ],
    );
  }

  async function favorite(poem: Poem) {
    await toggleFavorite(poem.id, !poem.favorite);
    await refresh();
  }

  function runClusters() {
    setClusters(clusterPoems(corpus, k));
  }

  const poemById = useMemo(
    () => new Map(corpus.map((poem) => [poem.id, poem])),
    [corpus],
  );

  return (
    <SafeAreaView style={styles.safe}>
      <StatusBar style="light" backgroundColor={win95.blue} />
      <View style={styles.desktop}>
        <Window title="PoemSpace 95 — ארכיון שירה">
          <View style={styles.menuBar}>
            <Text style={styles.menuItem}>קובץ</Text>
            <Text style={styles.menuItem}>עריכה</Text>
            <Text style={styles.menuItem}>תצוגה</Text>
            <Text style={styles.menuItem}>עזרה</Text>
          </View>

          <View style={styles.nav}>
            <RetroButton active={screen === "archive"} onPress={() => setScreen("archive")}>
              ארכיון
            </RetroButton>
            <RetroButton active={screen === "editor"} onPress={newPoem}>
              כתיבה
            </RetroButton>
            <RetroButton active={screen === "clusters"} onPress={() => setScreen("clusters")}>
              קלאסטרים
            </RetroButton>
            <RetroButton active={screen === "stats"} onPress={() => setScreen("stats")}>
              נתונים
            </RetroButton>
          </View>

          <View style={styles.content}>
            {screen === "archive" && (
              <View style={styles.fill}>
                <View style={styles.toolbar}>
                  <View style={styles.grow}>
                    <RetroField
                      value={query}
                      onChangeText={setQuery}
                      placeholder="חיפוש בכותרת, בטקסט או בתגיות"
                    />
                  </View>
                  <RetroButton onPress={newPoem}>+ חדש</RetroButton>
                </View>

                <View style={styles.toolbar}>
                  <View style={styles.monthField}>
                    <RetroField
                      value={month}
                      onChangeText={setMonth}
                      placeholder="חודש: YYYY-MM"
                    />
                  </View>
                  <RetroButton
                    active={sort === "written_desc"}
                    onPress={() => setSort("written_desc")}
                  >
                    חדש←ישן
                  </RetroButton>
                  <RetroButton
                    active={sort === "written_asc"}
                    onPress={() => setSort("written_asc")}
                  >
                    ישן←חדש
                  </RetroButton>
                </View>

                {error ? <Text style={styles.error}>{error}</Text> : null}

                <FlatList
                  data={poems}
                  keyExtractor={(item) => String(item.id)}
                  contentContainerStyle={poems.length ? styles.list : styles.emptyList}
                  ListEmptyComponent={
                    <View style={styles.empty}>
                      <Text style={styles.emptyTitle}>הארכיון ריק</Text>
                      <Text style={styles.rtlText}>לחץ "+ חדש" כדי להכניס את השיר הראשון.</Text>
                    </View>
                  }
                  renderItem={({ item }) => (
                    <Pressable onPress={() => editPoem(item)} style={styles.poemRow}>
                      <View style={styles.poemMain}>
                        <Text style={styles.poemTitle}>
                          {item.favorite ? "★ " : ""}
                          {item.title}
                        </Text>
                        <Text numberOfLines={2} style={styles.preview}>
                          {item.body}
                        </Text>
                        <Text style={styles.meta}>
                          {item.writtenAt ?? "ללא תאריך"}
                          {item.tags ? `   |   ${item.tags}` : ""}
                        </Text>
                      </View>
                      <View style={styles.rowActions}>
                        <RetroButton onPress={() => void favorite(item)}>
                          {item.favorite ? "☆" : "★"}
                        </RetroButton>
                        <RetroButton onPress={() => confirmDelete(item)}>מחק</RetroButton>
                      </View>
                    </Pressable>
                  )}
                />
              </View>
            )}

            {screen === "editor" && (
              <KeyboardAvoidingView
                behavior={Platform.OS === "ios" ? "padding" : undefined}
                style={styles.fill}
              >
                <ScrollView contentContainerStyle={styles.editor}>
                  <Text style={styles.label}>כותרת</Text>
                  <RetroField
                    value={draft.title}
                    onChangeText={(title) => setDraft((d) => ({ ...d, title }))}
                    placeholder="ללא כותרת"
                  />

                  <Text style={styles.label}>תאריך כתיבה</Text>
                  <RetroField
                    value={draft.writtenAt ?? ""}
                    onChangeText={(writtenAt) =>
                      setDraft((d) => ({ ...d, writtenAt: writtenAt || null }))
                    }
                    placeholder="YYYY-MM-DD"
                  />

                  <Text style={styles.label}>תגיות</Text>
                  <RetroField
                    value={draft.tags}
                    onChangeText={(tags) => setDraft((d) => ({ ...d, tags }))}
                    placeholder="למשל: לילה, בית, זיכרון"
                  />

                  <Text style={styles.label}>השיר</Text>
                  <RetroField
                    value={draft.body}
                    onChangeText={(body) => setDraft((d) => ({ ...d, body }))}
                    placeholder="כתוב או הדבק כאן..."
                    multiline
                  />

                  <View style={styles.editorButtons}>
                    <RetroButton
                      onPress={() => setDraft((d) => ({ ...d, favorite: !d.favorite }))}
                      active={draft.favorite}
                    >
                      ★ מועדף
                    </RetroButton>
                    <RetroButton onPress={() => void persistDraft()} disabled={busy}>
                      {busy ? "שומר..." : "שמור"}
                    </RetroButton>
                    <RetroButton onPress={() => setScreen("archive")}>ביטול</RetroButton>
                  </View>
                </ScrollView>
              </KeyboardAvoidingView>
            )}

            {screen === "clusters" && (
              <ScrollView contentContainerStyle={styles.analysis}>
                <Text style={styles.sectionTitle}>קלאסטרים מקומיים</Text>
                <Text style={styles.rtlText}>
                  TF-IDF + cosine k-means. כל החישוב מתבצע על המכשיר.
                </Text>

                <View style={styles.kRow}>
                  {[2, 3, 4, 5, 6].map((value) => (
                    <RetroButton key={value} active={k === value} onPress={() => setK(value)}>
                      K={value}
                    </RetroButton>
                  ))}
                </View>

                <RetroButton onPress={runClusters} disabled={!corpus.length}>
                  הרץ קלאסטרים על {corpus.length} שירים
                </RetroButton>

                {clusters.map((cluster) => (
                  <View key={cluster.id} style={styles.clusterBox}>
                    <View style={styles.clusterTitleBar}>
                      <Text style={styles.clusterTitle}>
                        אשכול {cluster.id} — {cluster.poemIds.length} שירים
                      </Text>
                    </View>
                    <Text style={styles.rtlText}>
                      מילים מובילות: {cluster.topTerms.join(" · ") || "—"}
                    </Text>
                    {cluster.poemIds.map((id) => {
                      const poem = poemById.get(id);
                      if (!poem) return null;
                      return (
                        <Pressable key={id} onPress={() => editPoem(poem)}>
                          <Text style={styles.poemLink}>• {poem.title}</Text>
                        </Pressable>
                      );
                    })}
                  </View>
                ))}
              </ScrollView>
            )}

            {screen === "stats" && (
              <ScrollView contentContainerStyle={styles.analysis}>
                <Text style={styles.sectionTitle}>סטטיסטיקות קורפוס</Text>

                <View style={styles.statsGrid}>
                  <Stat label="שירים" value={stats.poems} />
                  <Stat label="מילים" value={stats.words} />
                  <Stat label="מילים ייחודיות" value={stats.uniqueWords} />
                  <Stat label="ממוצע מילים לשיר" value={stats.averageWords} />
                </View>

                <View style={styles.infoBox}>
                  <Text style={styles.rtlText}>שיר מתוארך ראשון: {stats.firstDate ?? "—"}</Text>
                  <Text style={styles.rtlText}>שיר מתוארך אחרון: {stats.lastDate ?? "—"}</Text>
                </View>

                <View style={styles.infoBox}>
                  <Text style={styles.label}>מילים שכיחות בקורפוס</Text>
                  {stats.topTerms.map((item) => (
                    <View key={item.term} style={styles.termRow}>
                      <Text style={styles.termCount}>{item.count}</Text>
                      <Text style={styles.termName}>{item.term}</Text>
                    </View>
                  ))}
                </View>
              </ScrollView>
            )}
          </View>

          <View style={styles.statusBar}>
            <Text style={styles.statusText}>PoemSpace 95</Text>
            <Text style={styles.statusText}>
              {corpus.length} שירים | Offline
            </Text>
          </View>
        </Window>
      </View>
    </SafeAreaView>
  );
}

function Stat({ label, value }: { label: string; value: number }) {
  return (
    <View style={styles.statBox}>
      <Text style={styles.statValue}>{value}</Text>
      <Text style={styles.statLabel}>{label}</Text>
    </View>
  );
}

const styles = StyleSheet.create({
  safe: { flex: 1, backgroundColor: win95.desktop },
  desktop: { flex: 1, backgroundColor: win95.desktop, padding: 6 },
  window: {
    flex: 1,
    backgroundColor: win95.face,
    borderTopWidth: 2,
    borderLeftWidth: 2,
    borderTopColor: win95.light,
    borderLeftColor: win95.light,
    borderBottomWidth: 2,
    borderRightWidth: 2,
    borderBottomColor: win95.dark,
    borderRightColor: win95.dark,
  },
  titleBar: {
    height: 30,
    backgroundColor: win95.blue,
    flexDirection: "row-reverse",
    alignItems: "center",
    justifyContent: "space-between",
    paddingHorizontal: 4,
  },
  titleBarText: { color: win95.blueText, fontWeight: "700", fontSize: 14 },
  titleButtons: { flexDirection: "row", gap: 3 },
  titleMini: {
    width: 20,
    height: 20,
    textAlign: "center",
    color: win95.black,
    backgroundColor: win95.face,
    borderWidth: 1,
    borderTopColor: win95.light,
    borderLeftColor: win95.light,
    borderBottomColor: win95.dark,
    borderRightColor: win95.dark,
    fontWeight: "700",
  },
  windowBody: { flex: 1 },
  menuBar: {
    height: 28,
    flexDirection: "row-reverse",
    alignItems: "center",
    gap: 20,
    paddingHorizontal: 10,
    borderBottomWidth: 1,
    borderBottomColor: win95.mid,
  },
  menuItem: { color: win95.black, fontSize: 13 },
  nav: {
    flexDirection: "row-reverse",
    flexWrap: "wrap",
    gap: 5,
    padding: 6,
    borderBottomWidth: 2,
    borderBottomColor: win95.light,
  },
  button: {
    minHeight: 32,
    justifyContent: "center",
    alignItems: "center",
    paddingHorizontal: 10,
    backgroundColor: win95.face,
    borderTopWidth: 2,
    borderLeftWidth: 2,
    borderTopColor: win95.light,
    borderLeftColor: win95.light,
    borderBottomWidth: 2,
    borderRightWidth: 2,
    borderBottomColor: win95.dark,
    borderRightColor: win95.dark,
  },
  buttonPressed: {
    borderTopColor: win95.dark,
    borderLeftColor: win95.dark,
    borderBottomColor: win95.light,
    borderRightColor: win95.light,
  },
  buttonActive: {
    backgroundColor: win95.blue,
    borderTopColor: win95.dark,
    borderLeftColor: win95.dark,
    borderBottomColor: win95.light,
    borderRightColor: win95.light,
  },
  buttonText: { color: win95.black, fontWeight: "600", fontSize: 13 },
  buttonTextActive: { color: win95.blueText },
  disabled: { opacity: 0.5 },
  content: { flex: 1, minHeight: 0 },
  fill: { flex: 1 },
  toolbar: {
    flexDirection: "row-reverse",
    gap: 6,
    paddingHorizontal: 6,
    paddingTop: 6,
    alignItems: "center",
  },
  grow: { flex: 1 },
  monthField: { flex: 1, maxWidth: 180 },
  field: {
    minHeight: 38,
    backgroundColor: win95.field,
    color: win95.fieldText,
    paddingHorizontal: 8,
    paddingVertical: 7,
    fontSize: 15,
    writingDirection: "rtl",
    borderTopWidth: 2,
    borderLeftWidth: 2,
    borderTopColor: win95.dark,
    borderLeftColor: win95.dark,
    borderBottomWidth: 2,
    borderRightWidth: 2,
    borderBottomColor: win95.light,
    borderRightColor: win95.light,
  },
  textArea: {
    minHeight: 300,
    textAlignVertical: "top",
    lineHeight: 25,
    fontSize: 17,
  },
  list: { padding: 6, gap: 6 },
  emptyList: { flexGrow: 1, justifyContent: "center", padding: 20 },
  empty: {
    backgroundColor: win95.face,
    padding: 24,
    borderWidth: 2,
    borderTopColor: win95.dark,
    borderLeftColor: win95.dark,
    borderBottomColor: win95.light,
    borderRightColor: win95.light,
  },
  emptyTitle: {
    textAlign: "right",
    writingDirection: "rtl",
    fontSize: 18,
    fontWeight: "700",
    marginBottom: 8,
  },
  poemRow: {
    flexDirection: "row-reverse",
    alignItems: "stretch",
    gap: 8,
    padding: 8,
    backgroundColor: win95.face,
    borderTopWidth: 2,
    borderLeftWidth: 2,
    borderTopColor: win95.light,
    borderLeftColor: win95.light,
    borderBottomWidth: 2,
    borderRightWidth: 2,
    borderBottomColor: win95.dark,
    borderRightColor: win95.dark,
  },
  poemMain: { flex: 1 },
  poemTitle: {
    textAlign: "right",
    writingDirection: "rtl",
    color: win95.blue,
    fontSize: 17,
    fontWeight: "700",
  },
  preview: {
    textAlign: "right",
    writingDirection: "rtl",
    color: win95.black,
    lineHeight: 20,
    marginTop: 4,
  },
  meta: {
    textAlign: "right",
    writingDirection: "rtl",
    color: win95.dark,
    fontSize: 12,
    marginTop: 5,
  },
  rowActions: { gap: 5, justifyContent: "center" },
  editor: { padding: 10, gap: 6 },
  editorButtons: {
    flexDirection: "row-reverse",
    flexWrap: "wrap",
    gap: 8,
    marginTop: 6,
  },
  label: {
    textAlign: "right",
    writingDirection: "rtl",
    fontWeight: "700",
    color: win95.black,
    marginTop: 5,
  },
  analysis: { padding: 10, gap: 10 },
  sectionTitle: {
    textAlign: "right",
    writingDirection: "rtl",
    fontWeight: "700",
    fontSize: 20,
    color: win95.blue,
  },
  rtlText: {
    textAlign: "right",
    writingDirection: "rtl",
    color: win95.black,
    lineHeight: 21,
  },
  kRow: { flexDirection: "row-reverse", flexWrap: "wrap", gap: 5 },
  clusterBox: {
    backgroundColor: win95.face,
    borderWidth: 2,
    borderTopColor: win95.light,
    borderLeftColor: win95.light,
    borderBottomColor: win95.dark,
    borderRightColor: win95.dark,
    paddingBottom: 8,
  },
  clusterTitleBar: {
    backgroundColor: win95.blue,
    paddingHorizontal: 6,
    paddingVertical: 4,
    marginBottom: 7,
  },
  clusterTitle: {
    color: win95.blueText,
    textAlign: "right",
    writingDirection: "rtl",
    fontWeight: "700",
  },
  poemLink: {
    textAlign: "right",
    writingDirection: "rtl",
    color: win95.blue,
    paddingHorizontal: 8,
    paddingVertical: 3,
    textDecorationLine: "underline",
  },
  statsGrid: {
    flexDirection: "row-reverse",
    flexWrap: "wrap",
    gap: 8,
  },
  statBox: {
    minWidth: "47%",
    flexGrow: 1,
    backgroundColor: win95.face,
    padding: 12,
    borderWidth: 2,
    borderTopColor: win95.light,
    borderLeftColor: win95.light,
    borderBottomColor: win95.dark,
    borderRightColor: win95.dark,
  },
  statValue: {
    textAlign: "center",
    color: win95.blue,
    fontSize: 26,
    fontWeight: "700",
  },
  statLabel: { textAlign: "center", color: win95.black, marginTop: 4 },
  infoBox: {
    backgroundColor: win95.face,
    padding: 10,
    gap: 4,
    borderWidth: 2,
    borderTopColor: win95.dark,
    borderLeftColor: win95.dark,
    borderBottomColor: win95.light,
    borderRightColor: win95.light,
  },
  termRow: {
    flexDirection: "row-reverse",
    justifyContent: "space-between",
    borderBottomWidth: StyleSheet.hairlineWidth,
    borderBottomColor: win95.mid,
    paddingVertical: 4,
  },
  termName: { color: win95.black, writingDirection: "rtl" },
  termCount: { color: win95.dark },
  statusBar: {
    minHeight: 28,
    flexDirection: "row-reverse",
    justifyContent: "space-between",
    alignItems: "center",
    paddingHorizontal: 7,
    borderTopWidth: 2,
    borderTopColor: win95.light,
  },
  statusText: { color: win95.black, fontSize: 12 },
  error: {
    margin: 6,
    backgroundColor: "#ffffff",
    color: "#800000",
    padding: 8,
    textAlign: "right",
    writingDirection: "rtl",
  },
});
