package ge.dchechelashvili.notesapp.data

import android.content.Context
import androidx.room.Database
import androidx.room.Room
import androidx.room.RoomDatabase
import ge.dchechelashvili.notesapp.data.dao.ItemsDao
import ge.dchechelashvili.notesapp.data.entity.Item

@Database(entities = arrayOf(Item::class), version = 1)
abstract class ItemsDatabase(): RoomDatabase() {
    abstract fun itemsDao(): ItemsDao

    companion object{
        private val dbName = "items-db"
        private lateinit var INSTANCE: ItemsDatabase

        fun getInstance(): ItemsDatabase{
            return INSTANCE
        }

        fun createDatabase(context: Context){
            INSTANCE = Room.databaseBuilder(context, ItemsDatabase::class.java, dbName).build()
        }
    }
}