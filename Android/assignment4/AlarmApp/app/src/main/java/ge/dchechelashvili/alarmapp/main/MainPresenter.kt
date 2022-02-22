package ge.dchechelashvili.alarmapp.main

import android.app.AlarmManager
import android.app.PendingIntent
import android.content.Intent
import androidx.appcompat.app.AppCompatActivity
import ge.dchechelashvili.alarmapp.data.Alarm
import java.util.*

class MainPresenter(private var view: IMainView?) : IMainPresenter {
    private val interactor = MainInteractor(this)

    fun getAlarms(isFirstCall: Boolean = false) {
        interactor.getAlarmList(view, isFirstCall)
    }

    fun saveAlarms() {
        interactor.putAlarmList(view)
    }

    override fun onAlarmListFetched(alarms: List<Alarm>) {
        view?.showAlarmList(alarms)
    }

    override fun onModeFetched(mode: Int) {
        view?.setMode(mode)
    }

    fun removeAlarmItem(alarm: Alarm) {
        interactor.removeAlarmFromList(alarm)
    }

    fun addAlarm(hourOfDay: Int, minute: Int, view: AppCompatActivity) {
        val calendar = Calendar.getInstance()
        calendar.set(
            calendar.get(Calendar.YEAR), calendar.get(Calendar.MONTH), calendar.get(Calendar.DATE),
            hourOfDay, minute, 0
        )

        val time = calendar.timeInMillis - System.currentTimeMillis()

        if (time < 0) {
            return
        }

        var hour = hourOfDay.toString()
        var min = minute.toString()
        if (hourOfDay < 10) {
            hour = "0$hour"
        }
        if (minute < 10) {
            min = "0$min"
        }

        val alarm = Alarm("$hour:$min", true)
        interactor.addAlarmToList(alarm)

        activateAlarm(alarm, view, calendar)
    }


    fun getMode() {
        interactor.getAppMode(view)
    }

    fun saveMode(mode: Int) {
        interactor.putMode(view, mode)
    }

    fun detachView() {
        view = null
    }

    private fun activateAlarm(alarm: Alarm, view: AppCompatActivity, calendar: Calendar) {
        val (alarmManager, pi) = getAlarmManager(alarm, view)
        alarmManager.setExact(AlarmManager.RTC_WAKEUP, calendar.timeInMillis, pi)
    }

    fun deactivateAlarm(alarm: Alarm, view: AppCompatActivity) {
        val (alarmManager, pi) = getAlarmManager(alarm, view)
        alarmManager.cancel(pi)
    }

    private fun getAlarmManager(
        alarm: Alarm,
        view: AppCompatActivity
    ): Pair<AlarmManager, PendingIntent> {
        val requestCode =
            alarm.time.substringBefore(":").toInt() * 60 + alarm.time.substringAfter(":").toInt()
        val pi = PendingIntent.getBroadcast(
            view,
            requestCode,
            Intent(AlarmReceiver.ALARM_ACTION_NAME).apply {
                `package` = view.packageName
                putExtra("time", alarm.time)
            },
            0
        )
        val alarmManager = view.getSystemService(AppCompatActivity.ALARM_SERVICE) as AlarmManager
        return Pair(alarmManager, pi)
    }

    fun changeAlarmStatus(alarm: Alarm, checked: Boolean, view: AppCompatActivity) {
        if (checked) {
            val hourOfDay = alarm.time.substringBefore(":").toInt()
            val minute = alarm.time.substringAfter(":").toInt()
            val calendar = Calendar.getInstance()
            calendar.set(
                calendar.get(Calendar.YEAR),
                calendar.get(Calendar.MONTH),
                calendar.get(Calendar.DATE),
                hourOfDay,
                minute,
                0
            )
            activateAlarm(alarm, view, calendar)
        } else {
            deactivateAlarm(alarm, view)
        }
        interactor.updateAlarm(alarm, checked)
        saveAlarms()
    }

}