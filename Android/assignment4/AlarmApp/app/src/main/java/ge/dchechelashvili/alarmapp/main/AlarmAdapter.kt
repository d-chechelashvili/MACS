package ge.dchechelashvili.alarmapp.main

import android.util.Log
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import android.widget.TextView
import androidx.appcompat.widget.SwitchCompat
import androidx.recyclerview.widget.RecyclerView
import ge.dchechelashvili.alarmapp.R
import ge.dchechelashvili.alarmapp.data.Alarm

class AlarmAdapter(private val listListener: AlarmListListener) :
    RecyclerView.Adapter<AlarmAdapter.AlarmViewHolder>() {

    var list = listOf<Alarm>()

    override fun onCreateViewHolder(parent: ViewGroup, viewType: Int): AlarmViewHolder {
        return AlarmViewHolder(
            LayoutInflater.from(parent.context).inflate(R.layout.alarm_item, parent, false)
        )
    }

    override fun onBindViewHolder(holder: AlarmViewHolder, position: Int) {
        val alarm = list[position]
        holder.bindAlarm(alarm)
        holder.switch.setOnCheckedChangeListener { _, isChecked ->
            listListener.onAlarmStatusChanged(alarm, isChecked)
        }
        holder.itemView.setOnLongClickListener {
            listListener.onAlarmItemClicked(alarm)
            true
        }
    }

    override fun getItemCount(): Int {
        return list.size
    }

    inner class AlarmViewHolder(view: View) : RecyclerView.ViewHolder(view) {
        fun bindAlarm(alarm: Alarm) {
            time.text = alarm.time
            switch.setOnCheckedChangeListener(null)
            switch.isChecked = alarm.active
        }

        private val time: TextView = view.findViewById(R.id.time)
        val switch: SwitchCompat = view.findViewById(R.id.alarmStatus)
    }
}