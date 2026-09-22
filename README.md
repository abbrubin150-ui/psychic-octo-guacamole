# PoemSpace 95

Offline-first Android poetry archive with a Windows 95-inspired UI.

## What works

- Hebrew RTL writing and reading.
- Local SQLite archive; no server is required.
- Create, edit, delete, favorite and tag poems.
- Search by title, text or tags.
- Filter by month (`YYYY-MM`) and sort by writing date.
- Local TF-IDF + cosine k-means clustering.
- Corpus statistics and frequent terms.
- GitHub Actions builds an installable Android APK automatically.

## Local development

Requirements: Node.js 22.13+.

```bash
npm install
npx expo start
```

## APK

Every push to `main` runs **Build Android APK**.

In GitHub:

1. Open **Actions**.
2. Open the latest **Build Android APK** run.
3. Download the artifact named `poemspace95-android-apk`.
4. Extract `app-debug.apk` and install it on Android.

This is a debug-signed APK intended for direct installation/testing. No Expo account or cloud token is required for this build path.

## Architecture

```text
App.tsx
 ├─ src/db.ts       SQLite persistence
 ├─ src/analysis.ts TF-IDF, k-means, statistics
 ├─ src/types.ts    domain types
 └─ src/theme.ts    Windows 95 palette
```

The analysis engine and database are deliberately independent of the UI so they can be tested or replaced later without rewriting the app.
