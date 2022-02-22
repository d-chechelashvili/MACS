package ge.dchechelashvili.weatherapp

import androidx.fragment.app.Fragment
import androidx.fragment.app.FragmentActivity
import androidx.viewpager2.adapter.FragmentStateAdapter
import ge.dchechelashvili.weatherapp.fragments.HourlyForecastFragment
import ge.dchechelashvili.weatherapp.fragments.TodayForecastFragment

class ViewPagerAdapter(activity: FragmentActivity, private val fragmentsList: ArrayList<Fragment>): FragmentStateAdapter(activity) {

    override fun getItemCount(): Int {
       return fragmentsList.size
    }

    override fun createFragment(position: Int): Fragment {
        return fragmentsList[position]
    }

    fun updateTodayWeather(data: WeatherResponseModel){
        (fragmentsList[0] as TodayForecastFragment).update(data)
    }

    fun updateHourlyWeather(data: ForecastResponseModel){
        (fragmentsList[1] as HourlyForecastFragment).update(data)
    }
}