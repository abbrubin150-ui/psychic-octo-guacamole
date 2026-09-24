package com.pixelphysics.sandbox;

import android.app.Activity;
import android.os.Bundle;
import android.view.Window;
import android.view.WindowManager;

public class AndroidLauncher extends Activity {
    private PixelPhysicsVoxelView gameView;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        requestWindowFeature(Window.FEATURE_NO_TITLE);
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_FULLSCREEN);
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);

        gameView = new PixelPhysics25DView(this);
        setContentView(gameView);
    }

    @Override
    protected void onPause() {
        super.onPause();
        if (gameView != null) gameView.setKeepScreenOn(false);
    }

    @Override
    protected void onResume() {
        super.onResume();
        if (gameView != null) {
            gameView.setKeepScreenOn(true);
            gameView.requestFocus();
        }
    }
}
