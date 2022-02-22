package ge.dchechelashvili.weatherapp

import okhttp3.RequestBody
import retrofit2.Call
import retrofit2.http.GET
import retrofit2.http.Query

interface OpenWeatherMap {

    @GET("weather")
    fun getTodayWeather(@Query("q") cityName: String, @Query("appid") apiKey: String, @Query("units") units: String): Call<WeatherResponseModel>


    @GET("forecast")
    fun getForecast(@Query("q") cityName: String, @Query("appid") apiKey: String, @Query("units") units: String): Call<ForecastResponseModel>
}