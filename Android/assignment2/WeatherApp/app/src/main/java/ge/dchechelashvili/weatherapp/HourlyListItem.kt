package ge.dchechelashvili.weatherapp

import android.widget.ImageView

data class HourlyListItem(var dateLabel: String, var icon: String,
                      var temperature: String, var forecastLabel: String)