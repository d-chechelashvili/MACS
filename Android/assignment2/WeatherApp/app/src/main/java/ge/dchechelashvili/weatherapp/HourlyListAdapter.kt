package ge.dchechelashvili.weatherapp

import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import android.widget.ImageView
import android.widget.TextView
import androidx.recyclerview.widget.RecyclerView
import com.bumptech.glide.Glide

class HourlyListAdapter(var list: List<HourlyListItem>): RecyclerView.Adapter<HourlyItemViewHolder>() {

    override fun onCreateViewHolder(parent: ViewGroup, viewType: Int): HourlyItemViewHolder {
        val view = LayoutInflater.from(parent.context).inflate(R.layout.hourly_list_item, parent, false)
        return HourlyItemViewHolder(view)
    }

    override fun onBindViewHolder(holder: HourlyItemViewHolder, position: Int) {
        holder.dateLabel.text = list[position].dateLabel
        val icon = list[position].icon
        val imageUrl = "https://openweathermap.org/img/wn/$icon@2x.png"
        Glide.with(holder.itemView).load(imageUrl).into(holder.weatherImage)
        holder.temperature.text = list[position].temperature
        holder.forecastLabel.text = list[position].forecastLabel
    }

    override fun getItemCount(): Int {
        return list.size
    }

    fun updateData(data: List<HourlyListItem>){
        list = data
        notifyDataSetChanged()
    }
}

class HourlyItemViewHolder(itemView: View): RecyclerView.ViewHolder(itemView){
    var dateLabel = itemView.findViewById<TextView>(R.id.dateLabel)
    var weatherImage = itemView.findViewById<ImageView>(R.id.weatherImage)
    var temperature = itemView.findViewById<TextView>(R.id.item_temperature)
    var forecastLabel = itemView.findViewById<TextView>(R.id.item_description)
}