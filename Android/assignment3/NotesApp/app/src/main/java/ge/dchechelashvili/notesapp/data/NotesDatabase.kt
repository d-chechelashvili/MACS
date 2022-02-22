package ge.dchechelashvili.notesapp.data

import android.content.Context
import androidx.room.Database
import androidx.room.Room
import androidx.room.RoomDatabase
import ge.dchechelashvili.notesapp.data.dao.NotesDao
import ge.dchechelashvili.notesapp.data.entity.Note

@Database(entities = arrayOf(Note::class), version = 1)
abstract class NotesDatabase(): RoomDatabase() {
    abstract fun notesDao(): NotesDao

    companion object{
        private val dbName = "notes-db"
        private lateinit var INSTANCE: NotesDatabase

        fun getInstance(): NotesDatabase{
            return INSTANCE
        }

        fun createDatabase(context: Context){
            INSTANCE = Room.databaseBuilder(context, NotesDatabase::class.java, dbName).build()
        }
    }
}