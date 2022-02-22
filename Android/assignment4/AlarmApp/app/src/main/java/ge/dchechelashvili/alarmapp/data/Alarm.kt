package ge.dchechelashvili.alarmapp.data

data class Alarm(val time: String, val active: Boolean): Comparable<Alarm> {
    override fun compareTo(other: Alarm): Int {
        if (time > other.time) return 1
        if (time < other.time) return -1
        if (active && !other.active) return 1
        if (!active && other.active) return -1
        return 0
    }
}