package com.pixelphysics.sandbox;

import android.content.Context;
import android.content.SharedPreferences;
import android.os.Build;
import android.os.VibrationEffect;
import android.os.Vibrator;
import android.util.Log;
import android.webkit.JavascriptInterface;

import org.json.JSONObject;

public final class AndroidBridge {
    private static final String TAG = "PixelPhysicsV8";
    private static final int MAX_SAVE_BYTES = 2_000_000;

    private final Context context;
    private final SharedPreferences saves;
    private final Vibrator vibrator;

    AndroidBridge(Context context) {
        this.context = context.getApplicationContext();
        this.saves = this.context.getSharedPreferences("pixel_physics_enterprise_v8", Context.MODE_PRIVATE);
        this.vibrator = (Vibrator) this.context.getSystemService(Context.VIBRATOR_SERVICE);
    }

    @JavascriptInterface
    public void vibrate(int durationMs, int amplitude) {
        if (vibrator == null || !vibrator.hasVibrator()) return;
        int d = Math.max(1, Math.min(durationMs, 250));
        int a = Math.max(1, Math.min(amplitude, 255));
        try {
            if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.O) {
                vibrator.vibrate(VibrationEffect.createOneShot(d, a));
            } else {
                vibrator.vibrate(d);
            }
        } catch (Throwable t) {
            Log.w(TAG, "vibration failed", t);
        }
    }

    @JavascriptInterface
    public boolean saveState(String slot, String json) {
        if (!validSlot(slot) || json == null || json.length() > MAX_SAVE_BYTES) return false;
        try {
            // Validate before persisting so corrupt state never replaces a good save.
            new JSONObject(json);
            saves.edit()
                    .putString("save." + slot, json)
                    .putLong("save." + slot + ".ts", System.currentTimeMillis())
                    .apply();
            return true;
        } catch (Throwable t) {
            Log.e(TAG, "saveState rejected", t);
            return false;
        }
    }

    @JavascriptInterface
    public String loadState(String slot) {
        if (!validSlot(slot)) return "";
        return saves.getString("save." + slot, "");
    }

    @JavascriptInterface
    public void clearState(String slot) {
        if (!validSlot(slot)) return;
        saves.edit().remove("save." + slot).remove("save." + slot + ".ts").apply();
    }

    @JavascriptInterface
    public void log(String level, String message) {
        String m = message == null ? "" : message;
        switch (level == null ? "" : level) {
            case "error": Log.e(TAG, m); break;
            case "warn": Log.w(TAG, m); break;
            default: Log.i(TAG, m); break;
        }
    }

    @JavascriptInterface
    public void reportMetric(String name, double value) {
        if (name == null || name.length() > 80 || !Double.isFinite(value)) return;
        Log.i(TAG, "metric " + name + "=" + value);
    }

    @JavascriptInterface
    public String buildInfo() {
        try {
            JSONObject o = new JSONObject();
            o.put("version", "0.8.1-voxel-dynamics");
            o.put("sdk", Build.VERSION.SDK_INT);
            o.put("device", Build.DEVICE);
            o.put("model", Build.MODEL);
            return o.toString();
        } catch (Throwable t) {
            return "{}";
        }
    }

    private static boolean validSlot(String slot) {
        return slot != null && slot.matches("[A-Za-z0-9_.-]{1,48}");
    }
}
