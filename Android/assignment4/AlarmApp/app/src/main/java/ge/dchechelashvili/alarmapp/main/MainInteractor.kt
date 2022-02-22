package ge.dchechelashvili.alarmapp.main

import android.content.Context
import androidx.appcompat.app.AppCompatActivity
import ge.dchechelashvili.alarmapp.data.Alarm
import java.util.*

class MainInteractor(private val presenter: IMainPresenter) {

    private lateinit var alarms: SortedSet<Alarm>

    companion object {
        private const val PREF_NAME = "alarm-pref"
        private const val ALARM_LIST = "alarmList"
        private const val APP_MODE = "appMode"
    }

    fun getAlarmList(view: IMainView?, isFirstCall: Boolean = false) {
        if (isFirstCall) {
            val sharedPreferences =
                (view as AppCompatActivity).getSharedPreferences(PREF_NAME, Context.MODE_PRIVATE)
            val alarms = sharedPreferences.getStringSet(ALARM_LIST, mutableSetOf<String>())
            val alarmList = mutableSetOf<Alarm>()
            if (alarms != null) {
                for (alarm in alarms) {
                    alarmList.add(Alarm(alarm.dropLast(4), alarm.takeLast(4).toBoolean()))
                }
            }
            this.alarms = alarmList.toSortedSet()
        }
        presenter.onAlarmListFetched(alarms.toList())
    }

    fun putAlarmList(view: IMainView?) {
        val sharedPreferences =
            (view as AppCompatActivity).getSharedPreferences(PREF_NAME, Context.MODE_PRIVATE)
        val set = mutableSetOf<String>()
        for (alarm in alarms) {
            if (alarm.active) {
                set.add(alarm.time + "true")
            } else {
                set.add(alarm.time + "fals")
            }
        }
        sharedPreferences.edit().putStringSet(ALARM_LIST, set).apply()
    }

    fun getAppMode(view: IMainView?) {
        val sharedPreferences =
            (view as AppCompatActivity).getSharedPreferences(PREF_NAME, Context.MODE_PRIVATE)
        val mode = sharedPreferences.getInt(APP_MODE, -1)
        presenter.onModeFetched(mode)
    }

    fun putMode(view: IMainView?, mode: Int) {
        val sharedPreferences =
            (view as AppCompatActivity).getSharedPreferences(PREF_NAME, Context.MODE_PRIVATE)
        sharedPreferences.edit().putInt(APP_MODE, mode).apply()
    }

    fun removeAlarmFromList(alarm: Alarm) {
        alarms.remove(alarm)
        presenter.onAlarmListFetched(alarms.toList())
    }

    fun addAlarmToList(alarm: Alarm) {
        alarms.add(alarm)
        presenter.onAlarmListFetched(alarms.toList())
    }

    fun updateAlarm(alarm: Alarm, checked: Boolean) {
        alarms.add(Alarm(alarm.time, checked))
        removeAlarmFromList(alarm)
    }
}