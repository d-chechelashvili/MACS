package ge.dchechelashvili.notesapp.main

import android.util.Log
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import android.widget.CheckBox
import android.widget.TextView
import androidx.recyclerview.widget.RecyclerView
import ge.dchechelashvili.notesapp.R
import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

class NoteAdapter(
    val listListener: NoteListListener,
    var noteList: List<Note>,
    var items: List<List<Item>>
) :
    RecyclerView.Adapter<NoteAdapter.NoteViewHolder>() {

    override fun onCreateViewHolder(parent: ViewGroup, viewType: Int): NoteAdapter.NoteViewHolder {
        return NoteViewHolder(
            LayoutInflater.from(parent.context).inflate(R.layout.note_item, parent, false)
        )
    }

    override fun onBindViewHolder(holder: NoteAdapter.NoteViewHolder, position: Int) {
        val note = noteList[position]
        val noteItems = items[position]
        holder.bindNote(note, noteItems)
        holder.itemView.setOnClickListener {
            listListener.onNoteItemClicked(note)
        }
    }

    fun updateData(notes: List<Note>, items: List<List<Item>>) {
        noteList = notes
        this.items = items
        notifyDataSetChanged()
    }

    override fun getItemCount(): Int {
        return noteList.size
    }

    inner class NoteViewHolder(view: View) : RecyclerView.ViewHolder(view) {
        fun bindNote(note: Note, noteItems: List<Item>) {
            header.text = note.name

            for (i in 0..2) {
                items[i].visibility = View.GONE
                if (i >= noteItems.size)
                    continue
                val item = noteItems[i]
                items[i].isChecked = item.checked
                items[i].text = item.content
                items[i].visibility = View.VISIBLE
            }

            if (noteItems.size > 3) {
                more.visibility = View.VISIBLE
                itemCount.visibility = View.VISIBLE
                val str = "+ ${noteItems.size} checked items"
                itemCount.text = str
            } else {
                more.visibility = View.GONE
                itemCount.visibility = View.GONE
            }
        }

        private val header = view.findViewById<TextView>(R.id.header)
        private val items = listOf<CheckBox>(
            view.findViewById(R.id.item1),
            view.findViewById(R.id.item2), view.findViewById(R.id.item3)
        )
        private val more = view.findViewById<TextView>(R.id.more)
        private val itemCount = view.findViewById<TextView>(R.id.itemCount)
    }

}