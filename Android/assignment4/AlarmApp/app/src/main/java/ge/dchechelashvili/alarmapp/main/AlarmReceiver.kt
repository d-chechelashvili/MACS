package ge.dchechelashvili.alarmapp.main

import android.app.NotificationChannel
import android.app.NotificationManager.IMPORTANCE_HIGH
import android.app.PendingIntent
import android.content.BroadcastReceiver
import android.content.Context
import android.content.Intent
import android.os.Handler
import android.os.Looper
import androidx.core.app.NotificationCompat
import androidx.core.app.NotificationManagerCompat
import ge.dchechelashvili.alarmapp.R

class AlarmReceiver: BroadcastReceiver() {
    override fun onReceive(context: Context?, intent: Intent?) {
        context?.let {
            val notificationManager = NotificationManagerCompat.from(context)
            val time = intent?.getStringExtra("time")
            val requestCode =
                (time?.substringBefore(":")?.toInt() ?: 0) * 60 +
                        (time?.substringAfter(":")?.toInt() ?: 0)

            if (intent != null) {
                if (intent.action == CANCEL_ACTION_NAME) {
                    notificationManager.cancel(requestCode)
                    return
                }
            }

            val notificationClickPendingIntent = PendingIntent.getActivity(
                context,
                requestCode,
                Intent(context, MainActivity::class.java),
                0
            )

            val cancelButtonClick = PendingIntent.getBroadcast(
                context,
                requestCode,
                Intent(CANCEL_ACTION_NAME).apply {
                    `package` = context.packageName
                    putExtra("time", time)
                },
                0
            )

            val snoozeButtonClick = PendingIntent.getBroadcast(
                context,
                requestCode,
                Intent(SNOOZE_ACTION_NAME).apply {
                    `package` = context.packageName
                    putExtra("time", time)
                },
                0
            )

            createChannel(notificationManager)
            val notification = NotificationCompat.Builder(context, CHANNEL_ID)
                .setSmallIcon(R.drawable.alarm_clock)
                .setContentTitle("Alarm Message!")
                .setContentText("Alarm set on $time")
                .setContentIntent(notificationClickPendingIntent)
                .addAction(R.mipmap.ic_launcher, "CANCEL", cancelButtonClick)
                .addAction(R.mipmap.ic_launcher, "SNOOZE", snoozeButtonClick)
                .build()


            var delayMillis: Long = 0
            if (intent?.action == SNOOZE_ACTION_NAME) {
                delayMillis = 60000
                notificationManager.cancel(requestCode)
            }

            Handler(Looper.getMainLooper()).postDelayed({
                notificationManager.notify(requestCode, notification)
            }, delayMillis)
        }
    }

    private fun createChannel(notificationManager: NotificationManagerCompat) {
        val notificationChannel = NotificationChannel(CHANNEL_ID, "channel_name", IMPORTANCE_HIGH)
        notificationManager.createNotificationChannel(notificationChannel)
    }

    companion object{
        const val CANCEL_ACTION_NAME = "ge.dchechelashvili.alarmapp.main.CANCEL_ACTION"
        const val SNOOZE_ACTION_NAME = "ge.dchechelashvili.alarmapp.main.SNOOZE_ACTION"
        const val ALARM_ACTION_NAME = "ge.dchechelashvili.alarmapp.main.ALARM_ACTION"
        const val CHANNEL_ID = "ge.dchechelashvili.alarmapp.main.CHANNEL_1"
    }
}