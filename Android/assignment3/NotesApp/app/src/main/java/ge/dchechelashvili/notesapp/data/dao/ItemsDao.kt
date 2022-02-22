package ge.dchechelashvili.notesapp.data.dao

import androidx.room.*
import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

@Dao
interface ItemsDao {
    @Insert
    fun addItem(item: Item)

    @Insert(onConflict = OnConflictStrategy.REPLACE)
    fun updateItem(item: Item)

    @Query("SELECT * FROM item WHERE noteId = :noteId")
    fun getNoteItems(noteId: Int): List<Item>

    @Query("DELETE FROM item WHERE noteId = :noteId")
    fun removeItems(noteId: Int)
}