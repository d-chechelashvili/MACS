package ge.dchechelashvili.notesapp.data.entity

import androidx.room.ColumnInfo
import androidx.room.Entity
import androidx.room.PrimaryKey

@Entity
data class Item(
    @ColumnInfo val content: String,
    @ColumnInfo val noteId: Int,
    @ColumnInfo val checked: Boolean
){
    @PrimaryKey(autoGenerate = true)
    var id: Int = 0
}