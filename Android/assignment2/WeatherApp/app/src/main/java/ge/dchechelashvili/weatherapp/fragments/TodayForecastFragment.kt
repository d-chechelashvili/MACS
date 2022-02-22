package ge.dchechelashvili.weatherapp.fragments

import android.os.Bundle
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import android.widget.ImageView
import android.widget.TextView
import androidx.fragment.app.Fragment
import com.bumptech.glide.Glide
import ge.dchechelashvili.weatherapp.ForecastModel
import ge.dchechelashvili.weatherapp.R
import ge.dchechelashvili.weatherapp.WeatherResponseModel
import java.util.*
import kotlin.math.roundToInt

class TodayForecastFragment : Fragment() {

    lateinit var weatherIcon: ImageView
    lateinit var temperature: TextView
    lateinit var description: TextView
    lateinit var temperature2: TextView
    lateinit var feelsLike: TextView
    lateinit var humidity: TextView
    lateinit var pressure: TextView


    override fun onCreateView( inflater: LayoutInflater, container: ViewGroup?,
                               savedInstanceState: Bundle? ): View? {
        val view = inflater.inflate(R.layout.fragment_today_forecast, container, false)
        weatherIcon = view.findViewById<ImageView>(R.id.weatherIcon)
        temperature = view.findViewById<TextView>(R.id.temperature)
        description = view.findViewById<TextView>(R.id.description)
        temperature2 = view.findViewById<TextView>(R.id.temperature2)
        feelsLike = view.findViewById<TextView>(R.id.feelsLike)
        humidity = view.findViewById<TextView>(R.id.humidity)
        pressure = view.findViewById<TextView>(R.id.pressure)
        return view
    }

    fun update(data: WeatherResponseModel){
        val (timezone, main, weather) = data
        val (description, icon) = weather[0]
        val (temp, feels_like, pressure, humidity) = main
        val imageUrl = "https://openweathermap.org/img/wn/$icon@2x.png"
        Glide.with(this).load(imageUrl).into(weatherIcon)
        temperature.text  = convertToText(temp, "°")
        this.description.text = description.toUpperCase(Locale.ENGLISH)
        temperature2.text = convertToText(temp, "°")
        feelsLike.text = convertToText(feels_like, "°")
        this.humidity.text = convertToText(humidity, "%")
        this.pressure.text = pressure.toInt().toString()
    }

    private fun convertToText(toBeText: Double, addition: String): String{
        return toBeText.toInt().toString() + addition
    }
}