package ge.dchechelashvili.weatherapp.fragments

import android.os.Bundle
import android.util.Log
import androidx.fragment.app.Fragment
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import androidx.core.content.ContextCompat
import androidx.core.content.ContextCompat.getColor
import androidx.core.content.ContextCompat.getDrawable
import androidx.recyclerview.widget.DividerItemDecoration
import androidx.recyclerview.widget.RecyclerView
import ge.dchechelashvili.weatherapp.*
import java.text.DateFormatSymbols
import java.text.SimpleDateFormat
import java.time.format.DateTimeFormatter
import java.util.*

class HourlyForecastFragment : Fragment() {

    lateinit var rvHourly: RecyclerView

    override fun onCreateView(
        inflater: LayoutInflater, container: ViewGroup?,
        savedInstanceState: Bundle?
    ): View? {
        val view = inflater.inflate(R.layout.fragment_hourly_forecast, container, false)
        rvHourly = view.findViewById(R.id.rv_hourly)
        rvHourly.adapter = HourlyListAdapter(emptyList())
        val decoration = DividerItemDecoration(context, DividerItemDecoration.VERTICAL)
        getDrawable(view.context, R.drawable.divider_white)?.let { decoration.setDrawable(it) }
        rvHourly.addItemDecoration(decoration)
        return view
    }

    fun update(data: ForecastResponseModel) {
        val items = mutableListOf<HourlyListItem>()
        for (hourlyForecast in data.list){
            val date = adjustDate(hourlyForecast.dt_txt)
            val icon = hourlyForecast.weather[0].icon
            val temperature = hourlyForecast.main.temp.toInt().toString() + "°"
            val description = hourlyForecast.weather[0].description
            val item = HourlyListItem(date, icon, temperature, description)
            items.add(item)
        }
        (rvHourly.adapter as HourlyListAdapter).updateData(items)
    }

    private fun adjustDate(date: String): String {
        var newDate = ""
        var hour = date.substring(11, 13).toInt()
        if (hour > 12){
            hour -= 12
            newDate += "$hour PM "
        }else{
            newDate += "$hour AM "
        }
        val day = date.substring(8,10)
        newDate += "$day "
        val month = DateFormatSymbols().months[date.substring(5,7).toInt()-1].substring(0,3)
        newDate += month

        if (!newDate[1].isDigit()){
            newDate = "0$newDate"
        }
        if(!newDate[7].isDigit()){
            newDate.substring(0,6) + "0" + newDate.substring(6, newDate.length)
        }
        return newDate
    }

}