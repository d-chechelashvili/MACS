package ge.dchechelashvili.weatherapp

import android.text.format.DateUtils
import java.util.*

data class WeatherModel(val description: String, val icon: String)

data class MainModel(val temp: Double, val feels_like: Double, val pressure: Double, val humidity: Double)

data class ForecastResponseModel(val list: List<ForecastModel>)

data class WeatherResponseModel(val timezone: Double, val main: MainModel, val weather: List<WeatherModel>)

data class ForecastModel(val dt_txt: String, val main: MainModel, val weather: List<WeatherModel>)
