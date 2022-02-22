package ge.dchechelashvili.alarmapp.main

import android.app.TimePickerDialog
import android.os.Bundle
import android.view.View
import android.widget.TextView
import android.widget.TimePicker
import androidx.appcompat.app.AlertDialog
import androidx.appcompat.app.AppCompatActivity
import androidx.appcompat.app.AppCompatDelegate
import androidx.core.app.NotificationManagerCompat
import androidx.recyclerview.widget.RecyclerView
import ge.dchechelashvili.alarmapp.R
import ge.dchechelashvili.alarmapp.data.Alarm
import java.util.*

class MainActivity : AppCompatActivity(), IMainView, AlarmListListener,
    TimePickerDialog.OnTimeSetListener {
    private lateinit var rvAlarms: RecyclerView
    private lateinit var presenter: MainPresenter
    private var adapter = AlarmAdapter(this)

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)

        presenter = MainPresenter(this)
        presenter.getMode()
        initView()
        presenter.getAlarms(true)
    }

    override fun onDestroy() {
        presenter.detachView()
        super.onDestroy()
    }

    private fun initView() {
        rvAlarms = findViewById(R.id.rvAlarms)
        rvAlarms.adapter = adapter

        val modeSwitcher = findViewById<TextView>(R.id.modeSwitcher)
        if (AppCompatDelegate.getDefaultNightMode() == AppCompatDelegate.MODE_NIGHT_YES) {
            modeSwitcher.text = getString(R.string.switch_to_light)
        } else {
            modeSwitcher.text = getString(R.string.switch_to_dark)
        }
    }

    override fun showAlarmList(alarms: List<Alarm>) {
        adapter.list = alarms
        adapter.notifyDataSetChanged()
    }

    override fun setMode(mode: Int) {
        if (mode == 0) {
            AppCompatDelegate.setDefaultNightMode(AppCompatDelegate.MODE_NIGHT_NO)
        } else if (mode == 1) {
            AppCompatDelegate.setDefaultNightMode(AppCompatDelegate.MODE_NIGHT_YES)
        }
    }

    override fun onAlarmItemClicked(alarm: Alarm) {
        val dialog = AlertDialog.Builder(this)
            .setTitle("Are you sure you want to delete item?")
            .setPositiveButton("Yes"
            ) { dialog, _ ->
                presenter.removeAlarmItem(alarm)
                presenter.deactivateAlarm(alarm, this)
                presenter.saveAlarms()
                dialog?.dismiss()
            }
            .setNegativeButton("No") { dialog, _ -> dialog?.dismiss() }
            .create()
        dialog.show()
    }

    override fun onAlarmStatusChanged(alarm: Alarm, isChecked: Boolean) {
        presenter.changeAlarmStatus(alarm, isChecked, this)
    }

    fun addAlarm(view: View) {
        val calendar = Calendar.getInstance()
        val timePickerDialog = TimePickerDialog(
            this,
            this,
            calendar.get(Calendar.HOUR_OF_DAY),
            calendar.get(Calendar.MINUTE),
            true
        )
        timePickerDialog.show()
    }

    fun changeTheme(view: View) {
        if (AppCompatDelegate.getDefaultNightMode() == AppCompatDelegate.MODE_NIGHT_YES) {
            AppCompatDelegate.setDefaultNightMode(AppCompatDelegate.MODE_NIGHT_NO)
            presenter.saveMode(0)
        } else {
            AppCompatDelegate.setDefaultNightMode(AppCompatDelegate.MODE_NIGHT_YES)
            presenter.saveMode(1)
        }
    }

    override fun onTimeSet(view: TimePicker?, hourOfDay: Int, minute: Int) {
        presenter.addAlarm(hourOfDay, minute, this)
        presenter.saveAlarms()
    }
}

interface AlarmListListener {
    fun onAlarmItemClicked(alarm: Alarm)
    fun onAlarmStatusChanged(alarm: Alarm, isChecked: Boolean)
}