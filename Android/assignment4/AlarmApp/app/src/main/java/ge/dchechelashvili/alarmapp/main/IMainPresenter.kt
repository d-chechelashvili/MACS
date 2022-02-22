package ge.dchechelashvili.alarmapp.main

import ge.dchechelashvili.alarmapp.data.Alarm

interface IMainPresenter {
    fun onAlarmListFetched(alarms: List<Alarm>)
    fun onModeFetched(mode: Int)
}