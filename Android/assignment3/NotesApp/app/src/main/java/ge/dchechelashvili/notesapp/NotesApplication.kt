package ge.dchechelashvili.notesapp

import android.app.Application
import ge.dchechelashvili.notesapp.data.ItemsDatabase
import ge.dchechelashvili.notesapp.data.NotesDatabase

class NotesApplication: Application() {
    override fun onCreate() {
        super.onCreate()
        this.deleteDatabase("notes-db")
        this.deleteDatabase("items-db")
        ItemsDatabase.createDatabase(this)
        NotesDatabase.createDatabase(this)
    }
}