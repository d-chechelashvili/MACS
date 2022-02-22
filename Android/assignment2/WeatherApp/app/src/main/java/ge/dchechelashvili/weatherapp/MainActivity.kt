package ge.dchechelashvili.weatherapp

import android.os.Bundle
import android.util.Log
import android.view.View
import android.widget.TextView
import android.widget.Toast
import androidx.appcompat.app.AppCompatActivity
import androidx.viewpager2.widget.ViewPager2
import com.google.android.material.bottomnavigation.BottomNavigationView
import ge.dchechelashvili.weatherapp.fragments.HourlyForecastFragment
import ge.dchechelashvili.weatherapp.fragments.TodayForecastFragment
import retrofit2.Call
import retrofit2.Callback
import retrofit2.Response
import retrofit2.Retrofit
import retrofit2.converter.gson.GsonConverterFactory
import java.util.*

class MainActivity : AppCompatActivity() {

    private lateinit var viewPager: ViewPager2
    private var fragmentsList = arrayListOf(TodayForecastFragment(), HourlyForecastFragment())
    private lateinit var selectedCity: String

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)
        selectedCity = getString(R.string.tbilisi)
        setupFragments()
    }

    private fun setupFragments() {
        viewPager = findViewById(R.id.viewPager)
        viewPager.adapter = ViewPagerAdapter(this, fragmentsList)

        val navBar = findViewById<BottomNavigationView>(R.id.navBar)
        navBar.setOnNavigationItemSelectedListener {
            when (it.itemId){
                R.id.ic_today -> viewPager.currentItem = 0
                R.id.ic_hourly -> viewPager.currentItem = 1
            }
            true
        }

        viewPager.registerOnPageChangeCallback(object: ViewPager2.OnPageChangeCallback() {
            override fun onPageSelected(position: Int) {
                when (position){
                    0 -> navBar.selectedItemId = R.id.ic_today
                    1 -> navBar.selectedItemId = R.id.ic_hourly
                }
                updateWeather()
            }
        })
    }

    fun loadForecast(view: View) {
        setCityName(view.id)
        updateWeather()
    }

    private fun setCityName(id: Int) {
        when (id){
            R.id.GeorgiaImage-> selectedCity = getString(R.string.tbilisi)
            R.id.UKImage -> selectedCity = getString(R.string.london)
            R.id.JamaicaImage -> selectedCity = getString(R.string.kingston)
        }
        findViewById<TextView>(R.id.cityName).text = selectedCity
    }

    private fun updateWeather() {
        val retrofit = Retrofit.Builder().baseUrl(getString(R.string.weatherBaseUrl))
            .addConverterFactory(GsonConverterFactory.create()).build()
        val weatherMapApi = retrofit.create(OpenWeatherMap:: class.java)
        when(viewPager.currentItem){
            0 -> updateTodayWeather(weatherMapApi)
            1 -> updateHourlyWeather(weatherMapApi)
        }
    }

    private fun updateTodayWeather(weatherMapApi: OpenWeatherMap?) {
        val weatherCall = weatherMapApi?.getTodayWeather(selectedCity,
            getString(R.string.apiKey), getString(R.string.units))

        weatherCall?.enqueue(object: Callback<WeatherResponseModel>{
            override fun onResponse(call: Call<WeatherResponseModel>, response: Response<WeatherResponseModel>) {
                if (response.isSuccessful){
                    response.body()?.let {
                        changeBackgroundColor(it)
                        (viewPager.adapter as ViewPagerAdapter).updateTodayWeather(it)
                    }
                }else{
                    Toast.makeText(applicationContext, response.message(), Toast.LENGTH_SHORT).show()
                    //Log.d("error_today_weather", response.code().toString() + ": " +  response.message())
                }

            }

            override fun onFailure(call: Call<WeatherResponseModel>, t: Throwable) {
                Log.d("fail_today_weather", "FAILURE")
            }
        })
    }

    private fun changeBackgroundColor(forecast: WeatherResponseModel) {
        val calendar = Calendar.getInstance(TimeZone.getTimeZone("UTC"))
        val timezone = forecast.timezone.toLong() * 1000
        calendar.timeInMillis += timezone
        val time = calendar.get(Calendar.HOUR_OF_DAY)
        if (time in 7..18){
            window.decorView.setBackgroundColor(getColor(R.color.dayColor))
        }else{
            window.decorView.setBackgroundColor(getColor(R.color.nightColor))
        }
    }

    private fun updateHourlyWeather(weatherMapApi: OpenWeatherMap?) {
        val forecastCall = weatherMapApi?.getForecast(selectedCity,
            getString(R.string.apiKey), getString(R.string.units))

        forecastCall?.enqueue(object: Callback<ForecastResponseModel>{
            override fun onResponse(call: Call<ForecastResponseModel>, response: Response<ForecastResponseModel>) {
                if (response.isSuccessful){
                    response.body()?.let {
                        (viewPager.adapter as ViewPagerAdapter).updateHourlyWeather(it)
                    }
                }else{
                    Toast.makeText(applicationContext, response.message(), Toast.LENGTH_SHORT).show()
                    //Log.d("error_forecast", response.code().toString() + ": " +  response.message())
                }
            }

            override fun onFailure(call: Call<ForecastResponseModel>, t: Throwable) {
                Log.e("fail_forecast", "FAILURE")
            }
        })
    }
}