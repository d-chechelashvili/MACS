package ge.dchechelashvili.notesapp.main

import SpacesItemDecoration
import android.content.Context
import android.content.Intent
import androidx.appcompat.app.AppCompatActivity
import android.os.Bundle
import android.util.Log
import android.view.KeyEvent
import android.view.View
import android.widget.Button
import android.widget.EditText
import android.widget.TextView
import androidx.recyclerview.widget.RecyclerView
import ge.dchechelashvili.notesapp.R
import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note
import ge.dchechelashvili.notesapp.details.DetailsActivity

class MainActivity : AppCompatActivity(), IMainView, NoteListListener {

    private lateinit var presenter: MainPresenter
    private lateinit var rvNotesPinned: RecyclerView
    private lateinit var rvNotesUnpinned: RecyclerView
    private var pinnedAdapter = NoteAdapter(this, emptyList(), emptyList())
    private var unpinnedAdapter = NoteAdapter(this, emptyList(), emptyList())
    private lateinit var searchBar: EditText

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)
        initView()
        presenter = MainPresenter(this)
        presenter.getNotes("%%")
    }

    override fun onDestroy() {
        presenter.detachView()
        super.onDestroy()
    }

    private fun initView() {
        rvNotesPinned = findViewById(R.id.rvNotesPinned)
        rvNotesPinned.adapter = pinnedAdapter
        rvNotesUnpinned = findViewById(R.id.rvNotesUnpinned)
        rvNotesUnpinned.adapter = unpinnedAdapter
        val verticalSpacing = resources.getDimensionPixelSize(R.dimen.vertical_margin)
        val horizontalSpacing = resources.getDimensionPixelSize(R.dimen.horizontal_margin)
        rvNotesPinned.addItemDecoration(SpacesItemDecoration(verticalSpacing, horizontalSpacing))
        rvNotesUnpinned.addItemDecoration(SpacesItemDecoration(verticalSpacing, horizontalSpacing))

        searchBar = findViewById(R.id.editTextTextSearchNotes)
        searchBar.setOnKeyListener(View.OnKeyListener { v, keyCode, event ->
            if (keyCode == KeyEvent.KEYCODE_ENTER && event.action == KeyEvent.ACTION_UP) {
                presenter.getNotes("%${searchBar.text.toString()}%")
                return@OnKeyListener true
            }
            false
        })
    }


    override fun showPinnedNotes(pinnedNotes: List<Note>, pinnedItems: List<List<Item>>) {
        if (pinnedNotes.isEmpty()){
            findViewById<TextView>(R.id.pinned).visibility = View.GONE
            findViewById<TextView>(R.id.other).visibility = View.GONE
            rvNotesPinned.visibility = View.GONE
        }else{
            findViewById<TextView>(R.id.pinned).visibility = View.VISIBLE
            findViewById<TextView>(R.id.other).visibility = View.VISIBLE
            rvNotesPinned.visibility = View.VISIBLE
        }
        pinnedAdapter.updateData(pinnedNotes, pinnedItems)
    }

    override fun showUnpinnedNotes(unpinnedNotes: List<Note>, unpinnedItems: List<List<Item>>) {
        unpinnedAdapter.updateData(unpinnedNotes, unpinnedItems)
    }

    override fun onNoteItemClicked(note: Note) {
        startDetails(note.id)
    }

    fun startDetails(id: Int){
        val intent = Intent(Intent(this, DetailsActivity::class.java))
        intent.putExtra(DetailsActivity.NOTE_ID_EXTRA, id)
        startActivityForResult(intent, 1)
    }

    override fun onActivityResult(requestCode: Int, resultCode: Int, data: Intent?) {
        super.onActivityResult(requestCode, resultCode, data)
        if (requestCode==1){
            presenter.getNotes("%%")
        }
    }

    fun addNote(view: View) {
        startDetails(-1)
    }

}

interface NoteListListener {
    fun onNoteItemClicked(note: Note)
}