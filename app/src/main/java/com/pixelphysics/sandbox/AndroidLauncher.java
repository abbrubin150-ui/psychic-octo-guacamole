package com.pixelphysics.sandbox;

import android.os.Bundle;
import com.badlogic.gdx.backends.android.AndroidApplication;
import com.badlogic.gdx.backends.android.AndroidApplicationConfiguration;

public class AndroidLauncher extends AndroidApplication {
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        AndroidApplicationConfiguration cfg = new AndroidApplicationConfiguration();
        cfg.useImmersiveMode = false;
        cfg.useWakelock = false;
        cfg.useAccelerometer = false;
        cfg.useCompass = false;
        cfg.useGyroscope = false;
        cfg.disableAudio = true;

        initialize(new PixelPhysics25DGame(), cfg);
    }
}
