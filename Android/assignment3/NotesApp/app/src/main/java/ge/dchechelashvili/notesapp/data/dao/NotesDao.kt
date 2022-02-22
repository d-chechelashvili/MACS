package ge.dchechelashvili.notesapp.data.dao

import androidx.room.Dao
import androidx.room.Insert
import androidx.room.OnConflictStrategy
import androidx.room.Query
import ge.dchechelashvili.notesapp.data.entity.Note

@Dao
interface NotesDao {
    @Insert
    fun addNote(note: Note)

    @Insert(onConflict = OnConflictStrategy.REPLACE)
    fun updateNote(note: Note): Long

    @Query("SELECT * FROM Note WHERE name LIKE :filterString ORDER BY modificationDate")
    fun getNotes(filterString: String): List<Note>

    @Query("SELECT * FROM Note WHERE id=:id")
    fun getNote(id: Int): Note


    @Query("DELETE FROM Note WHERE id=:noteId")
    fun removeNote(noteId: Int)
}