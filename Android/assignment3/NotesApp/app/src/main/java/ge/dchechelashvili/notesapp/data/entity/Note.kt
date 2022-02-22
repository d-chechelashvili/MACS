package ge.dchechelashvili.notesapp.data.entity

import androidx.room.ColumnInfo
import androidx.room.Entity
import androidx.room.PrimaryKey

@Entity
data class Note(
    @ColumnInfo val name: String,
    @ColumnInfo val modificationDate: Long,
    @ColumnInfo val pinned: Boolean
){
    @PrimaryKey(autoGenerate = true)
    var id: Int = 0
}