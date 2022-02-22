package ge.dchechelashvili.alarmapp.main

import ge.dchechelashvili.alarmapp.data.Alarm

interface IMainView {
    fun showAlarmList(alarms: List<Alarm>)
    fun setMode(mode: Int)
}