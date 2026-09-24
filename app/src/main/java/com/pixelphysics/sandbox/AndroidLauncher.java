package com.pixelphysics.sandbox;

import android.app.Activity;
import android.graphics.Color;
import android.net.Uri;
import android.os.Bundle;
import android.util.Log;
import android.view.Gravity;
import android.view.Window;
import android.view.WindowManager;
import android.webkit.ConsoleMessage;
import android.webkit.RenderProcessGoneDetail;
import android.webkit.ConsoleMessage;
import android.webkit.WebChromeClient;
import android.webkit.WebResourceError;
import android.webkit.WebResourceRequest;
import android.webkit.WebResourceResponse;
import android.webkit.WebSettings;
import android.webkit.WebView;
import android.webkit.WebViewClient;
import android.widget.TextView;

import androidx.annotation.Nullable;
import androidx.webkit.WebViewAssetLoader;

public class AndroidLauncher extends Activity {
    private static final String TAG = "PixelPhysicsV8";
    private static final String APP_ORIGIN = "https://appassets.androidplatform.net";
    private WebView webView;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        requestWindowFeature(Window.FEATURE_NO_TITLE);
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_FULLSCREEN);
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);

        startWebRuntime();
    }

    private void startWebRuntime() {
        final WebViewAssetLoader assetLoader = new WebViewAssetLoader.Builder()
                .addPathHandler("/assets/", new WebViewAssetLoader.AssetsPathHandler(this))
                .build();

        webView = new WebView(this);
        webView.setBackgroundColor(Color.rgb(11, 14, 20));
        webView.setOverScrollMode(WebView.OVER_SCROLL_NEVER);

        WebSettings s = webView.getSettings();
        s.setJavaScriptEnabled(true);
        s.setDomStorageEnabled(true);
        s.setAllowFileAccess(false);
        s.setAllowContentAccess(false);
        s.setBuiltInZoomControls(false);
        s.setDisplayZoomControls(false);
        s.setSupportZoom(false);
        s.setMediaPlaybackRequiresUserGesture(true);
        s.setMixedContentMode(WebSettings.MIXED_CONTENT_NEVER_ALLOW);
        s.setCacheMode(WebSettings.LOAD_NO_CACHE);

        webView.addJavascriptInterface(new AndroidBridge(this), "AndroidBridge");
        webView.setWebChromeClient(new WebChromeClient() {
            @Override
            public boolean onConsoleMessage(ConsoleMessage message) {
                String text = "JS " + message.messageLevel() + " "
                        + message.message() + " @"
                        + message.sourceId() + ":" + message.lineNumber();
                switch (message.messageLevel()) {
                    case ERROR:
                        Log.e(TAG, text);
                        break;
                    case WARNING:
                        Log.w(TAG, text);
                        break;
                    default:
                        Log.i(TAG, text);
                        break;
                }
                return true;
            }
        });

        webView.setWebViewClient(new WebViewClient() {
            @Override
            public @Nullable WebResourceResponse shouldInterceptRequest(WebView view, WebResourceRequest request) {
                return assetLoader.shouldInterceptRequest(request.getUrl());
            }

            @Override
            public void onPageFinished(WebView view, String url) {
                android.util.Log.i("VoxelDynamics", "PAGE_FINISHED " + url);
                super.onPageFinished(view, url);
            }

            @Override
            public void onReceivedError(
                    WebView view,
                    WebResourceRequest request,
                    android.webkit.WebResourceError error) {
                android.util.Log.e(
                        "VoxelDynamics",
                        "WEB_ERROR " + request.getUrl() + " " + error.getDescription()
                );
                super.onReceivedError(view, request, error);
            }

            @Override
            public boolean shouldOverrideUrlLoading(WebView view, WebResourceRequest request) {
                Uri u = request.getUrl();
                return !("https".equals(u.getScheme())
                        && "appassets.androidplatform.net".equals(u.getHost()));
            }

            @Override
            public void onReceivedError(WebView view, WebResourceRequest request, WebResourceError error) {
                super.onReceivedError(view, request, error);
                Log.e(TAG, "WEB_ERROR code=" + error.getErrorCode()
                        + " desc=" + error.getDescription()
                        + " url=" + request.getUrl());
            }

            @Override
            public boolean onRenderProcessGone(WebView view, RenderProcessGoneDetail detail) {
                Log.e(TAG, "WEBVIEW_RENDERER_GONE crash=" + detail.didCrash());
                showNativeFatal("Renderer process stopped. Reopen Pixel Physics Lab.");
                return true;
            }
        });

        setContentView(webView);
        Log.i(TAG, "HOST_START Pixel Physics Lab 0.8.1");
        webView.loadUrl(APP_ORIGIN + "/assets/www/index.html");
    }

    private void showNativeFatal(String message) {
        if (webView != null) {
            try {
                webView.destroy();
            } catch (Throwable ignored) {
            }
            webView = null;
        }

        TextView fallback = new TextView(this);
        fallback.setBackgroundColor(Color.rgb(11, 14, 20));
        fallback.setTextColor(Color.rgb(255, 216, 210));
        fallback.setTextSize(18);
        fallback.setGravity(Gravity.CENTER);
        fallback.setPadding(40, 40, 40, 40);
        fallback.setText("PIXEL PHYSICS LAB\n\n" + message);
        setContentView(fallback);
    }

    @Override
    protected void onPause() {
        if (webView != null) {
            webView.evaluateJavascript("window.__onAppPause && window.__onAppPause()", null);
            webView.onPause();
        }
        super.onPause();
    }

    @Override
    protected void onResume() {
        super.onResume();
        if (webView != null) {
            webView.onResume();
            webView.evaluateJavascript("window.__onAppResume && window.__onAppResume()", null);
        }
    }

    @Override
    protected void onDestroy() {
        if (webView != null) {
            webView.loadUrl("about:blank");
            webView.removeAllViews();
            webView.destroy();
            webView = null;
        }
        super.onDestroy();
    }
}
